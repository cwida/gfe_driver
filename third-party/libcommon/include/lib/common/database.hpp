/**
 * This file is part of libcommon.
 *
 * libcommon is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, orF
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

#ifndef COMMON_DATABASE_HPP
#define COMMON_DATABASE_HPP

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "error.hpp"

namespace common {

// An error caused by an operation in the database
DEFINE_EXCEPTION(DatabaseError);


/**
 * A wrapper to SQLite3 to record the outcome of experiments. Sample usage:
 *
 * Database db("path/to/data.sqlite3");
 * db.create_execution()("algorithm", "btree").save();
 *
 * vector<pair<string, string>> params; // these will be stored in the table parameters
 * params.push_back(pair<string, string>{"block_size", "32"});
 * params.push_back(pair<string, string>{"leaf_size", "64"});
 * db.store_parameters(params);
 *
 * // store the results of an outcome
 * db.add("experiment_aging")("completion_time", 32);
 *
 *
 * The result will be a star with the following tables:
 * 1- executions: it records the execution id and acts as the fact/central table
 * 2- parameters: a list of key/value parameters associated to the underlying execution
 * 3- experiment_aging: a table with the final outcome results.
 */
class Database {
public:
    class Execution;

    friend class Execution; // forward declaration

private:
    Database(const Database&) = delete;

    Database& operator=(Database&) = delete;

    const std::string m_database_path; // the path to the database connection
    void* m_handle; // opaque handle, actual connection to the database
    using ExecutionPtr = std::shared_ptr<Execution>;
    std::vector<ExecutionPtr> m_executions; // list of executions still valid/running
    bool m_keep_alive; // whether to keep the connection opened to the database, after each operation

public: // Record
    enum FieldType { TYPE_TEXT, TYPE_INTEGER, TYPE_REAL };

    struct AbstractField {
        std::string key;
        FieldType type;

        AbstractField(const std::string& key, FieldType type);

        virtual ~AbstractField();
    };

    struct TextField : public AbstractField {
        std::string value;

        TextField(const std::string& key, const std::string& value);
    };

    struct IntegerField : public AbstractField {
        int64_t value;

        IntegerField(const std::string& key, int64_t value);
    };

    struct RealField : public AbstractField {
        double value;

        RealField(const std::string& key, double value);
    };

    class BaseRecord {
    protected:
        std::vector<std::shared_ptr<AbstractField>> m_fields;
    public:
        void add(const std::string& key, const std::string& value);
        void add(const std::string& key, int64_t value);
        void add(const std::string& key, uint64_t value);
        void add(const std::string& key, double value);
        void add(const BaseRecord& record);
        const std::vector<std::shared_ptr<AbstractField>>& fields() const;
        void dump(std::ostream& out) const;
    };

    template<typename Subclass>
    class Record : public BaseRecord {
    protected:
        void check_valid(){ /* nop */ }

    public:
        Subclass& operator()(const std::string& key, const std::string& value){
            reinterpret_cast<Subclass*>(this)->check_valid(); // callback to check the instance has not been already finalised
            add(key, value);
            return reinterpret_cast<Subclass&>(*this);
        }

        Subclass& operator()(const std::string& key, int64_t value){
            reinterpret_cast<Subclass*>(this)->check_valid();
            add(key, value);
            return reinterpret_cast<Subclass&>(*this);
        }

        Subclass& operator()(const std::string& key, double value){
            reinterpret_cast<Subclass*>(this)->check_valid();
            add(key, value);
            return reinterpret_cast<Subclass&>(*this);
        }

        // Treat all integer types as int64_t
        template<typename T>
        std::enable_if_t<std::is_integral_v<T>, Subclass&>
        operator()(const std::string& key, T value) {
            return operator()(key, static_cast<int64_t>(value));
        }

        template<typename Other>
        Subclass& operator()(const BaseRecord& other){
            reinterpret_cast<Subclass*>(this)->check_valid();
            add(other);
            return reinterpret_cast<Subclass&>(*this);
        }
    };


public: // Results
    friend class OutcomeBuilder;

    class OutcomeBuilder : public Record<OutcomeBuilder> {
        friend class Execution;
        OutcomeBuilder(OutcomeBuilder& object) = delete;
        OutcomeBuilder& operator=(OutcomeBuilder& object) = delete;

        std::shared_ptr<Execution> m_instance;
        const std::string m_table_name;

        Database* database() const;

        std::shared_ptr<Execution> execution() const;

        void save();
    public:
        OutcomeBuilder(std::shared_ptr<Execution> instance, const std::string& tableName);

        OutcomeBuilder(OutcomeBuilder&& object);

        ~OutcomeBuilder() noexcept(false);

        void dump(std::ostream& out) const;
    };

public: // Execution

    class ExecutionBuilder : public Record<ExecutionBuilder> {
        friend class Database;
        Database* m_instance;

        ExecutionBuilder(Database* instance);

    protected:
        void check_valid();

    public:
        ~ExecutionBuilder() noexcept(false);

        std::shared_ptr<Execution> save();
    };

    friend class Execution;

    class Execution {
        friend class ExecutionBuilder;

        Database* m_instance; // the wrapper associated to this instance
        const int64_t m_id; // SQLite3 id associated to this instance
        std::weak_ptr<Execution> m_self; // reference to the shared_ptr owning this instance

        Execution(Database* handle, int64_t execution_id);

    public:
        /**
         * Destructor. It automatically `closes' the current execution.
         */
        ~Execution();

        /**
         * Save into the table `parameters' a list of key/value arguments
         */
        void store_parameters(const std::vector<std::pair<std::string, std::string>>& params);

        /**
         * Add the given experiment results in the table `tableName'. The table
         * is created if it does not already exist.
         */
        OutcomeBuilder add(const std::string& tableName);

        OutcomeBuilder add(const char* tableName);

        /**
         * Retrieve the associated database
         */
        Database* database() const noexcept;

        /**
         * Retrieve the execution ID
         */
        int64_t id() const noexcept;

        /**
         * Check whether this execution is still valid
         */
        bool valid() const noexcept;

        /**
         * Terminate the current execution
         */
        void close();
    };


public:
    Database(const std::string& path, bool keep_alive = true);


    ~Database();

    /**
     * Create a new execution
     */
    ExecutionBuilder create_execution();

    /**
     * Open a connection to the SQLite3 database. This function becomes a `nop' if the wrapper is already connected
     */
    void connect();

    /**
     * Check whether the instance is currently connected to the database
     */
    bool is_connected() const noexcept;

    /**
    * Close the connection to the underlying SQLite3 database
    */
    void disconnect();

    /**
     * Retrieve the underlying connection handle
     */
    void* get_connection_handle() const noexcept;

    /**
     * The path to the database
     */
    const char* db_path() const noexcept;

    /**
     * Check whether the connection remains open after each operation
     */
    bool is_keep_alive() const noexcept;

    /**
     * Set the setting for keep alive
     * @param value true to maintain the connection open after each operation, false otherwise
     */
    void set_keep_alive(bool value) noexcept;

    /**
     * Retrieve the current execution.
     */
    std::shared_ptr<Execution> current() const;

    /**
     * Save into the table `parameters' a list of key/value arguments for the current execution
     */
    void store_parameters(const std::vector<std::pair<std::string, std::string>>& params);

    /**
     * Add the outcome results for the current execution in the table `table_name'. The table
     * is created if it does not already exist.
     */
    OutcomeBuilder add(const std::string& table_name);
};

} // namespace common

/**
 * Print the content of a record
 */
std::ostream& operator<<(std::ostream& out, const common::Database::BaseRecord& record);
std::ostream& operator<<(std::ostream& out, const common::Database::OutcomeBuilder& record);

#endif // COMMON_DATABASE_HPP
