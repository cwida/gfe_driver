/**
 * This file is part of libcommon.
 *
 * libcommon is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * libcommon is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with libcommon.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "database.hpp"

#include <cassert>
#include <exception> // std::uncaught_exceptions
#include <iostream>
#include <sqlite3.h>
#include <sstream>

#undef CURRENT_ERROR_TYPE
#define CURRENT_ERROR_TYPE DatabaseError

using namespace common;
using namespace std;

/*****************************************************************************
 *                                                                           *
 *   Debug                                                                   *
 *                                                                           *
 *****************************************************************************/
//#define DEBUG
#define COUT_DEBUG_FORCE(msg) cout << "[" << __FUNCTION__ << " @ " << __FILE__ << ":" << __LINE__ << "] " << msg << endl
#if defined(DEBUG)
#define COUT_DEBUG(msg) COUT_DEBUG_FORCE(msg)
#else
#define COUT_DEBUG(msg)
#endif

/*****************************************************************************
 *                                                                           *
 *  Database                                                                 *
 *                                                                           *
 *****************************************************************************/

Database::Database(const std::string &path, bool keep_alive) :
        m_database_path(path), m_handle(nullptr), m_keep_alive(keep_alive){
    // always attempt a connection on init
    connect();

    if(!is_keep_alive()) disconnect();
}

Database::~Database(){
    if(is_connected() || !m_executions.empty()){
        int rc {0};
        connect();
        auto connection = reinterpret_cast<sqlite3*>(m_handle);

        // Check there are no transactions active
        rc = sqlite3_get_autocommit(connection);
        if(rc == 0){ // Ignore the result, just attempt to rollback at this point
            sqlite3_exec(connection, "ROLLBACK", nullptr, nullptr, nullptr);
        }

        // Close all executions still active
        for(auto p_exec : m_executions){
            if(p_exec->valid())
                p_exec->close();
        }
        m_executions.clear();

        disconnect();
    }
}

void* Database::get_connection_handle() const noexcept {
    return m_handle;
}

void Database::connect(){
    if(is_connected()) return; // already connected?

    sqlite3* connection(nullptr);
    int rc = sqlite3_open(db_path(), &connection);
    if(rc != SQLITE_OK) { ERROR("Cannot open a SQLite connection to `" << m_database_path << "'"); }
    assert(connection != nullptr);
    m_handle = connection;
}

void Database::disconnect(){
    if(!is_connected()) return;

    int rc {0};
    auto connection = reinterpret_cast<sqlite3*>(m_handle);

    rc = sqlite3_close(connection);
    if(rc != SQLITE_OK){ // don't throw an exception here
        cerr << "[Database::disconnect] ERROR: " << sqlite3_errmsg(connection) << " error code: " << sqlite3_errcode(connection) << endl;
    }
    m_handle = nullptr;
}

bool Database::is_connected() const noexcept {
    return m_handle != nullptr;
}

bool Database::is_keep_alive() const noexcept {
    return m_keep_alive;
}

void Database::set_keep_alive(bool value) noexcept {
    m_keep_alive = value;
}

const char* Database::db_path() const noexcept {
    return m_database_path.c_str();
}

/*****************************************************************************
 *                                                                           *
 *  Connection                                                               *
 *                                                                           *
 *****************************************************************************/
namespace {
class Connection {
    Connection(const Connection&) = delete;
    Connection& operator=(Connection&) = delete;
    Database* m_instance;

public:
    Connection(Database* db) : m_instance(db){
        m_instance->connect();
    }

    ~Connection(){
        if(!m_instance->is_keep_alive())
            m_instance->disconnect();
    }

    operator sqlite3*(){
        void* handle = m_instance->get_connection_handle();
        assert(handle != nullptr);
        return reinterpret_cast<sqlite3*>(handle);
    }
};
} // anonymous namespace

/*****************************************************************************
 *                                                                           *
 *  Transaction                                                              *
 *                                                                           *
 *****************************************************************************/
 namespace {
 class Transaction{
     Connection* m_connection;

 public:
 Transaction(Connection& connection) : m_connection(nullptr) {
     // Start the transaction
     char* errmsg {nullptr};
     sqlite3* handle = connection;
     int rc = sqlite3_exec(handle, "BEGIN TRANSACTION", nullptr, nullptr, &errmsg);
     if(rc != SQLITE_OK || errmsg != nullptr){
         string error = errmsg; sqlite3_free(errmsg); errmsg = nullptr;
         ERROR("Cannot start the transaction: " << error);
     }

     m_connection = &connection;
 }

 ~Transaction(){
     if(m_connection != nullptr){
         if(uncaught_exceptions()){
             rollback();
         } else {
             commit();
         }
     }
 }

 void commit(){
     if(m_connection == nullptr) ERROR("Already closed");
     char* errmsg {nullptr};
     int rc = sqlite3_exec(*m_connection, "COMMIT", nullptr, nullptr, &errmsg);
     if(rc != SQLITE_OK || errmsg != nullptr){
         string error = errmsg; sqlite3_free(errmsg); errmsg = nullptr;
         ERROR("Cannot complete the transaction: " << error);
     }

     m_connection = nullptr;
 }


 void rollback(){
     if(m_connection == nullptr) ERROR("Already closed");
     sqlite3_exec(*m_connection, "ROLLBACK", nullptr, nullptr, nullptr);
     m_connection = nullptr;
 }

 };
 } // anonymous namespace



/*****************************************************************************
 *                                                                           *
 *  Record                                                                   *
 *                                                                           *
 *****************************************************************************/
Database::AbstractField::AbstractField(const string& key, FieldType type) : key(key), type(type) {
    // check the key is not id or exec_id
    locale loc;
    stringstream lstrb;
    for(auto e : key){ lstrb << tolower(e, loc); }
    string lstr = lstrb.str();
    if(lstr == "id" || lstr == "exec_id"){
        ERROR("Invalid attribute name: `" << key << "'. This name is reserved.");
    }
}

Database::AbstractField::~AbstractField() { }
Database::TextField::TextField(const string& key, const string& value) : AbstractField{ key, TYPE_TEXT }, value(value) { }
Database::IntegerField::IntegerField(const string& key, int64_t value) : AbstractField{ key, TYPE_INTEGER }, value(value) { }
Database::RealField::RealField(const string& key, double value) : AbstractField{ key, TYPE_REAL }, value(value) { }

void Database::BaseRecord::add(const string& key, const string& value) {
    shared_ptr<AbstractField> ptr(new TextField(key, value));
    m_fields.push_back(move(ptr));
}

void Database::BaseRecord::add(const string& key, int64_t value) {
    shared_ptr<AbstractField> ptr(new IntegerField(key, value));
    m_fields.push_back(move(ptr));
}

void Database::BaseRecord::add(const std::string& key, uint64_t value){
    add(key, static_cast<int64_t>(value));
}

void Database::BaseRecord::add(const string& key, double value) {
    shared_ptr<AbstractField> ptr(new RealField(key, value));
    m_fields.push_back(move(ptr));
}

void Database::BaseRecord::add(const BaseRecord& record) {
    for(auto ptr : record.fields()){
        m_fields.push_back(ptr);
    }
}

const vector<shared_ptr<Database::AbstractField>>& Database::BaseRecord::fields() const{
    return m_fields;
}

void Database::BaseRecord::dump(std::ostream& out) const{
    for(size_t i = 0; i < m_fields.size(); i++){
        auto e = m_fields[i];
        out << "[" << (i+1) << "] name: " << e->key << ", type: ";
        switch(e->type){
        case TYPE_TEXT:
            out << "text, value: \"" << dynamic_cast<TextField*>(e.get())->value << "\"";
            break;
        case TYPE_INTEGER:
            out << "int, value: " << dynamic_cast<IntegerField*>(e.get())->value;
            break;
        case TYPE_REAL:
            out << "real, value: " << dynamic_cast<RealField*>(e.get())->value;
            break;
        default:
            out << "unknown (" << e->type << ")";
        }
        out << "\n";
    }
}

/*****************************************************************************
 *                                                                           *
 *  Execution                                                                *
 *                                                                           *
 *****************************************************************************/
Database::ExecutionBuilder Database::create_execution(){
    return ExecutionBuilder{this};
}

shared_ptr<Database::Execution> Database::current() const {
    if(m_executions.empty()) ERROR("There are no executions registered");
    return m_executions.back();
}

void Database::store_parameters(const std::vector<std::pair<std::string, std::string>> &params) {
    current()->store_parameters(params);
}

Database::OutcomeBuilder Database::add(const std::string& table_name) {
    return current()->add(table_name);
}

Database::ExecutionBuilder::ExecutionBuilder(Database* instance) : m_instance(instance){ }

Database::ExecutionBuilder::~ExecutionBuilder() noexcept(false) {
    if(m_instance != nullptr) save();
}

void Database::ExecutionBuilder::check_valid() {
    if(m_instance == nullptr) ERROR("Instance already finalised");
}

std::shared_ptr<Database::Execution> Database::ExecutionBuilder::save() {
    check_valid();

    int rc = 0;
    Connection connection(m_instance);
    Transaction transaction(connection);
    const char* tableName = "executions";
    char* errmsg = nullptr;

    { // first check whether the table exists
        const char* SQL_check_table_exists = "SELECT 1 FROM sqlite_master WHERE type='table' AND name=?;";
        sqlite3_stmt* stmt (nullptr);
        rc = sqlite3_prepare_v2(connection, SQL_check_table_exists, -1, &stmt, nullptr);
        if(rc != SQLITE_OK || stmt == nullptr)
        ERROR("Cannot prepare the statement to insert an execution: " << sqlite3_errstr(rc));
        rc = sqlite3_bind_text(stmt, 1, tableName, /* compute length with strlen */ -1, /* do not free */ SQLITE_STATIC);
        if(rc != SQLITE_OK){ ERROR("Cannot bind the parameter: " << sqlite3_errstr(rc)); }
        rc = sqlite3_step(stmt);
        if(rc != SQLITE_DONE && rc != SQLITE_ROW){
            ERROR("Cannot execute the statement: " << sqlite3_errstr(rc));
        }
        bool table_exists = rc != SQLITE_DONE;
        rc = sqlite3_finalize(stmt); stmt = nullptr;
        assert(rc == SQLITE_OK);

        // the table `tableName' does not exist
        if(!table_exists){
            stringstream sqlcc;
            sqlcc << "CREATE TABLE " << tableName << "( ";
            sqlcc << "id INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT, ";
            sqlcc << "timeStart TIMESTAMP DEFAULT CURRENT_TIMESTAMP, ";
            sqlcc << "timeEnd TIMESTAMP ";
            for(auto& e : m_fields){
                sqlcc << ", " << e->key;
                switch(e->type){
                case TYPE_TEXT:
                    sqlcc << " TEXT NOT NULL"; break;
                case TYPE_INTEGER:
                    sqlcc << " INTEGER NOT NULL"; break;
                case TYPE_REAL:
                    sqlcc << " REAL NOT NULL"; break;
                default:
                ERROR("Invalid type: " << (int) e->type);
                }
            }
            sqlcc << ")";
            auto SQL_create_table = sqlcc.str();
            COUT_DEBUG("SQL DDL: " << SQL_create_table);
            rc = sqlite3_exec(connection, SQL_create_table.c_str(), nullptr, nullptr, &errmsg);
            if(rc != SQLITE_OK || errmsg != nullptr){
                string error = errmsg; sqlite3_free(errmsg); errmsg = nullptr;
                ERROR("Cannot create the table `" << tableName << "': " << error);
            }
        }
    } // does the table exist?

    { // Insert the results
        stringstream sqlcc;
        if(m_fields.size() > 0){
            sqlcc << "INSERT INTO " << tableName << " ( ";
            for(size_t i = 0; i < m_fields.size(); i++ ){
                if(i > 0) sqlcc << ", ";
                sqlcc << m_fields[i]->key;
            }
            sqlcc << " ) VALUES ( ";
            for(size_t i = 0; i < m_fields.size(); i++){
                if(i > 0) sqlcc << ", ";
                sqlcc << "?";
            }
            sqlcc << " )";
        } else { // m_fields.size() == 0
            sqlcc << "INSERT INTO " << tableName << " DEFAULT VALUES";
        }
        auto sqlccs = sqlcc.str();
        COUT_DEBUG("SQL DDL: " << sqlccs);

        sqlite3_stmt* stmt (nullptr);
        rc = sqlite3_prepare_v2(connection, sqlccs.c_str(), -1, &stmt, nullptr);
        if(rc != SQLITE_OK || stmt == nullptr){
            ERROR("Cannot prepare the statement to insert a result: " << sqlite3_errstr(rc));
        }
//        rc = sqlite3_bind_int64(stmt, 1, m_instance->id());
//        if(rc != SQLITE_OK){ ERROR("SQL Insert -> Results: cannot bind the parameter #1: " << sqlite3_errstr(rc)); }
        int index = 1;
        for(auto& e : m_fields){
            switch(e->type){
            case TYPE_TEXT:
                rc = sqlite3_bind_text(stmt, index, dynamic_cast<TextField*>(e.get())->value.c_str(), /* compute length with strlen */ -1, /* do not free */ SQLITE_STATIC);
                break;
            case TYPE_INTEGER:
                rc = sqlite3_bind_int64(stmt, index, dynamic_cast<IntegerField*>(e.get())->value);
                break;
            case TYPE_REAL:
                rc = sqlite3_bind_double(stmt, index, dynamic_cast<RealField*>(e.get())->value);
                break;
            default:
            ERROR("[Database::ResultsBuilder::save] Invalid type: " << (int) e->type);
            }
            if(rc != SQLITE_OK){ ERROR("SQL Insert -> Results: cannot bind the parameter " << index << ": " << sqlite3_errstr(rc)); }
            index++;
        }
        rc = sqlite3_step(stmt);
        if(rc != SQLITE_DONE){ ERROR("SQL Insert -> Results: cannot insert the values: " << sqlite3_errstr(rc) << ". SQL Statement: " << sqlccs); }
        rc = sqlite3_finalize(stmt); stmt = nullptr;
        assert(rc == SQLITE_OK);
    }

    // Retrieve the execution id
    auto execution_id = sqlite3_last_insert_rowid(connection);

    // Create the execution
    std::shared_ptr<Execution> execution;
    execution.reset(new Execution(m_instance, execution_id));
    execution->m_self = execution;
    m_instance->m_executions.push_back(execution);
    m_instance = nullptr; // avoid being called again!
    return execution;
}


Database::Execution::Execution(Database* handle, int64_t execution_id) : m_instance(handle), m_id(execution_id) {  }

Database::Execution::~Execution() {
    if(valid()) close();
}

int64_t Database::Execution::id() const noexcept {
    return m_id;
}

bool Database::Execution::valid() const noexcept {
    return m_instance != nullptr;
}

void Database::Execution::store_parameters(const vector<pair<string, string>>& params){
    if(!valid()) ERROR("Instance already terminated");

    Connection connection(m_instance);
    Transaction transaction(connection);
    int rc(0); char* errmsg {nullptr};

    // Create the table parameters, if it doesn't already exist
    auto SQL_create_table_parameters = ""
       "CREATE TABLE IF NOT EXISTS parameters ("
       "   exec_id INTEGER NOT NULL, "
       "   name TEXT NOT NULL, "
       "   value TEXT NOT NULL, "
       "   PRIMARY KEY(exec_id, name), "
       "   FOREIGN KEY(exec_id) REFERENCES executions ON DELETE CASCADE ON UPDATE CASCADE"
       ");";
    rc = sqlite3_exec(connection, SQL_create_table_parameters, nullptr, nullptr, &errmsg);
    if(rc != SQLITE_OK || errmsg != nullptr){
        string error = errmsg; sqlite3_free(errmsg); errmsg = nullptr;
        ERROR("Cannot create the table `parameters': " << error);
    }

    sqlite3_stmt* stmt { nullptr };
    auto SQL_insert_parameters = "INSERT INTO parameters (exec_id, name, value) VALUES (?, ?, ?);";
    rc = sqlite3_prepare_v2(connection, SQL_insert_parameters, -1, &stmt, nullptr);
    if(rc != SQLITE_OK || stmt == nullptr)
        ERROR("Cannot prepare the statement to insert the execution parameters: " << sqlite3_errstr(rc));

    for(auto p : params){
        const string& name = p.first;
        const string& value = p.second;

        rc = sqlite3_reset(stmt);
        if(rc != SQLITE_OK){ ERROR("SQL Insert -> Parameter [" << name << "]: cannot reset the statement: " << sqlite3_errstr(rc)); }
        rc = sqlite3_bind_int64(stmt, 1, id());
        if(rc != SQLITE_OK){ ERROR("SQL Insert -> Parameter [" << name << "]: cannot bind the execution ID: " << sqlite3_errstr(rc)); }
        rc = sqlite3_bind_text(stmt, 2, name.c_str(), -1, SQLITE_STATIC);
        if(rc != SQLITE_OK){ ERROR("SQL Insert -> Parameter [" << name << "]: cannot bind the key: " << sqlite3_errstr(rc)); }
        rc = sqlite3_bind_text(stmt, 3, value.c_str(), -1, SQLITE_STATIC);
        if(rc != SQLITE_OK){ ERROR("SQL Insert -> Parameter [" << name << "]: cannot bind the value: " << sqlite3_errstr(rc)); }
        rc = sqlite3_step(stmt);
        if(rc != SQLITE_DONE){ ERROR("SQL Insert -> Parameters: cannot insert the values for the parameter " << name << ": " << sqlite3_errstr(rc)); }
    }
    rc = sqlite3_finalize(stmt); stmt = nullptr;
    assert(rc == SQLITE_OK);
}

void Database::Execution::close(){
    if(!valid()) ERROR("Already terminated");
    Connection connection(m_instance);

    int rc (0);
    sqlite3_stmt* stmt (nullptr);
    auto SQL_update_execution = "UPDATE executions SET timeEnd = CURRENT_TIMESTAMP WHERE id = ?";



    rc = sqlite3_prepare_v2(connection, SQL_update_execution, -1, &stmt, nullptr);
    if(rc != SQLITE_OK || stmt == nullptr){
        cerr << "[Database::Execution::close] ERROR: Cannot prepare the statement to insert an execution: " << sqlite3_errstr(rc) << endl;
        goto next;
    }
    rc = sqlite3_bind_int64(stmt, 1, id());
    if(rc != SQLITE_OK){
        cerr << "[Database::Execution::close] ERROR: SQL Insert -> Executions: cannot bind the parameter #1: " << sqlite3_errstr(rc) << endl;
        goto next;
    }
    rc = sqlite3_step(stmt);
    if(rc != SQLITE_DONE){
        cerr << "[Database::Execution::close] ERROR: SQL Insert -> Executions: cannot insert the values: " << sqlite3_errstr(rc) << endl;
        goto next;
    }
    rc = sqlite3_finalize(stmt); stmt = nullptr;
    assert(rc == SQLITE_OK);

next:
    m_instance = nullptr;
}

Database* Database::Execution::database() const noexcept {
    return m_instance;
}

Database::OutcomeBuilder Database::Execution::add(const string& table_name) {
    return Database::OutcomeBuilder(m_self.lock(), table_name);
}

Database::OutcomeBuilder Database::Execution::add(const char* table_name) {
    string str_table_name{table_name};
    return add(str_table_name);
}

/*****************************************************************************
 *                                                                           *
 *  Outcomes                                                                 *
 *                                                                           *
 *****************************************************************************/
Database::OutcomeBuilder::OutcomeBuilder(std::shared_ptr<Execution> instance, const string& tableName) : m_instance(instance), m_table_name(tableName){
    if(instance.get() == nullptr || !instance->valid()){ ERROR("This execution has already been sealed"); }
}

Database::OutcomeBuilder::OutcomeBuilder(OutcomeBuilder&& object) : m_instance(object.m_instance), m_table_name(object.m_table_name){
    object.m_instance.reset();
}

Database::OutcomeBuilder::~OutcomeBuilder() noexcept(false) {
    if(m_instance) save();
}

shared_ptr<Database::Execution> Database::OutcomeBuilder::execution() const {
    if(m_instance.get() == nullptr || !m_instance->valid()) ERROR("Invalid execution instance");
    return m_instance;
}

Database* Database::OutcomeBuilder::database() const {
    auto e = execution();
    if(e->database() == nullptr) ERROR("Connection not available, the underlying execution has been invalidated");
    return e->database();
}

void Database::OutcomeBuilder::save(){ // this method can be invoked only by the dtor
    Connection connection(database());
    Transaction transaction(connection);
    int rc = 0;
    char* errmsg = nullptr;

    { // first check whether the table exists
        const char* SQL_check_table_exists = "SELECT 1 FROM sqlite_master WHERE type='table' AND name=?;";
        sqlite3_stmt* stmt (nullptr);
        rc = sqlite3_prepare_v2(connection, SQL_check_table_exists, -1, &stmt, nullptr);
        if(rc != SQLITE_OK || stmt == nullptr)
        ERROR("Cannot prepare the statement to insert an execution: " << sqlite3_errstr(rc));
        rc = sqlite3_bind_text(stmt, 1, m_table_name.c_str(), /* compute length with strlen */ -1, /* do not free */ SQLITE_STATIC);
        if(rc != SQLITE_OK){ ERROR("Cannot bind the parameter: " << sqlite3_errstr(rc)); }
        rc = sqlite3_step(stmt);
        if(rc != SQLITE_DONE && rc != SQLITE_ROW){
            ERROR("Cannot execute the statement: " << sqlite3_errstr(rc));
        }
        bool table_exists = rc != SQLITE_DONE;
        rc = sqlite3_finalize(stmt); stmt = nullptr;
        assert(rc == SQLITE_OK);

        // the table `tableName' does not exist
        if(!table_exists){
            stringstream sqlcc;
            sqlcc << "CREATE TABLE " << m_table_name << "( ";
            sqlcc << "id INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT, ";
            sqlcc << "exec_id INTEGER NOT NULL, ";
            for(auto e : m_fields){
                sqlcc << e->key;
                switch(e->type){
                case TYPE_TEXT:
                    sqlcc << " TEXT NOT NULL, "; break;
                case TYPE_INTEGER:
                    sqlcc << " INTEGER NOT NULL, "; break;
                case TYPE_REAL:
                    sqlcc << " REAL NOT NULL, "; break;
                default:
                    ERROR("Invalid type: " << (int) e->type);
                }
            }
            sqlcc << "FOREIGN KEY(exec_id) REFERENCES executions ON DELETE CASCADE ON UPDATE CASCADE";
            sqlcc << ")";
            auto SQL_create_table = sqlcc.str();
            rc = sqlite3_exec(connection, SQL_create_table.c_str(), nullptr, nullptr, &errmsg);
            if(rc != SQLITE_OK || errmsg != nullptr){
                string error = errmsg; sqlite3_free(errmsg); errmsg = nullptr;
                ERROR("Cannot create the table `" << m_table_name << "': " << error);
            }
        }
    } // does the table exist?

    { // Insert the results
        stringstream sqlcc;
        sqlcc << "INSERT INTO " << m_table_name << " ( exec_id";
        for(size_t i = 0; i < m_fields.size(); i++ ){
            sqlcc << ", " << m_fields[i]->key;
        }
        sqlcc << " ) VALUES ( ?";
        for(size_t i = 0; i < m_fields.size(); i++){
            sqlcc << ", ?";
        }
        sqlcc << " )";
        auto sqlccs = sqlcc.str();

        sqlite3_stmt* stmt (nullptr);
        rc = sqlite3_prepare_v2(connection, sqlccs.c_str(), -1, &stmt, nullptr);
        if(rc != SQLITE_OK || stmt == nullptr){
            ERROR("Cannot prepare the statement to insert a result: " << sqlite3_errstr(rc));
        }
        rc = sqlite3_bind_int64(stmt, 1, execution()->id());
        if(rc != SQLITE_OK){ ERROR("SQL Insert -> Results: cannot bind the exec_id: " << sqlite3_errstr(rc)); }
        int index = 2;
        for(auto& e : m_fields){
            switch(e->type){
            case TYPE_TEXT:
                rc = sqlite3_bind_text(stmt, index, dynamic_cast<TextField*>(e.get())->value.c_str(), /* compute length with strlen */ -1, /* do not free */ SQLITE_STATIC);
                break;
            case TYPE_INTEGER:
                rc = sqlite3_bind_int64(stmt, index, dynamic_cast<IntegerField*>(e.get())->value);
                break;
            case TYPE_REAL:
                rc = sqlite3_bind_double(stmt, index, dynamic_cast<RealField*>(e.get())->value);
                break;
            default:
            ERROR("[Database::ResultsBuilder::save] Invalid type: " << (int) e->type);
            }
            if(rc != SQLITE_OK){ ERROR("SQL Insert -> Results: cannot bind the parameter " << index << ": " << sqlite3_errstr(rc)); }
            index++;
        }
        rc = sqlite3_step(stmt);
        if(rc != SQLITE_DONE){ ERROR("SQL Insert -> Results: cannot insert the values: " << sqlite3_errstr(rc) << ". SQL Statement: " << sqlccs); }
        rc = sqlite3_finalize(stmt); stmt = nullptr;
        assert(rc == SQLITE_OK);
    }
}

void Database::OutcomeBuilder::dump(std::ostream& out) const{
    out << "table: " << m_table_name << ", # fields: " << fields().size() << "\n";
    BaseRecord::dump(out);
}

ostream& operator<<(ostream& out, const Database::BaseRecord& record) {
    record.dump(out);
    return out;
}

ostream& operator<<(ostream& out, const Database::OutcomeBuilder& instance) {
    instance.dump(out);
    return out;
}

