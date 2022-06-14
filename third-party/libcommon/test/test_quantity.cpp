#include "gtest/gtest.h"

#include "lib/common/quantity.hpp"

using namespace std;
using namespace common;

/**
 * Unit testing for the static method ComputerQuantity::parse()
 */
TEST(ComputerQuantity, parse){
    using ByteSuffix = ComputerQuantity::ByteSuffix;

    // straighforward numbers
    ASSERT_EQ(ComputerQuantity::parse("0", ByteSuffix::Missing), 0 );
    ASSERT_EQ(ComputerQuantity::parse("7", ByteSuffix::Missing), 7 );
    ASSERT_THROW(ComputerQuantity::parse("-1", ByteSuffix::Missing), QuantityError );

    // units
    ASSERT_EQ(ComputerQuantity::parse("0k", ByteSuffix::Missing), 0 );
    ASSERT_EQ(ComputerQuantity::parse("7k", ByteSuffix::Missing), 7 * 1024 );
    ASSERT_EQ(ComputerQuantity::parse("1K", ByteSuffix::Missing), 1024);
    ASSERT_EQ(ComputerQuantity::parse("1 k", ByteSuffix::Missing), 1024);
    ASSERT_THROW(ComputerQuantity::parse("-1k", ByteSuffix::Missing), QuantityError );
    ASSERT_EQ(ComputerQuantity::parse("1m", ByteSuffix::Missing), (1ull << 20));
    ASSERT_EQ(ComputerQuantity::parse("1M", ByteSuffix::Missing), (1ull << 20));
    ASSERT_EQ(ComputerQuantity::parse("1g", ByteSuffix::Missing), (1ull << 30));
    ASSERT_EQ(ComputerQuantity::parse("1G", ByteSuffix::Missing), (1ull << 30));
    ASSERT_EQ(ComputerQuantity::parse("1t", ByteSuffix::Missing), (1ull << 40));
    ASSERT_EQ(ComputerQuantity::parse("1T", ByteSuffix::Missing), (1ull << 40));
    ASSERT_THROW(ComputerQuantity::parse("0b", ByteSuffix::Missing), QuantityError );
    ASSERT_THROW(ComputerQuantity::parse("0B", ByteSuffix::Missing), QuantityError );
    ASSERT_THROW(ComputerQuantity::parse("0Kb", ByteSuffix::Missing), QuantityError );
    ASSERT_THROW(ComputerQuantity::parse("0 kB", ByteSuffix::Missing), QuantityError );
    ASSERT_THROW(ComputerQuantity::parse("0 Mb", ByteSuffix::Missing), QuantityError );
    ASSERT_THROW(ComputerQuantity::parse("0 MB", ByteSuffix::Missing), QuantityError );

    // byte quantities, with mandatory suffix
    ASSERT_THROW(ComputerQuantity::parse("0", ByteSuffix::Mandatory), QuantityError );
    ASSERT_THROW(ComputerQuantity::parse("0k", ByteSuffix::Mandatory), QuantityError );
    ASSERT_EQ(ComputerQuantity::parse("0b", ByteSuffix::Mandatory), 0 );
    ASSERT_EQ(ComputerQuantity::parse("0KB", ByteSuffix::Mandatory), 0 );
    ASSERT_EQ(ComputerQuantity::parse("2 KB", ByteSuffix::Mandatory), 2 * 1024 );

    // byte quantities, with optional byte suffix
    ASSERT_EQ(ComputerQuantity::parse("0", ByteSuffix::Optional), 0 );
    ASSERT_EQ(ComputerQuantity::parse("0k", ByteSuffix::Optional), 0 );
    ASSERT_EQ(ComputerQuantity::parse("0b", ByteSuffix::Optional), 0 );
    ASSERT_EQ(ComputerQuantity::parse("0KB", ByteSuffix::Optional), 0 );
    ASSERT_EQ(ComputerQuantity::parse("2 k", ByteSuffix::Optional), 2 * 1024 );
    ASSERT_EQ(ComputerQuantity::parse("2 KB", ByteSuffix::Optional), 2 * 1024 );
}

TEST(ComputerQuantity, to_string) {

    ASSERT_EQ(ComputerQuantity("1").to_string(), string("1"));
    ASSERT_EQ(ComputerQuantity("1 k").to_string(), string("1K"));
    ASSERT_EQ(ComputerQuantity("1 m ").to_string(), string("1M"));
    ASSERT_EQ(ComputerQuantity("1 G ").to_string(), string("1G"));
    ASSERT_EQ(ComputerQuantity("1", true).to_string(), string("1 bytes"));
    ASSERT_EQ(ComputerQuantity("1 k", true).to_string(), string("1 KB"));
    ASSERT_EQ(ComputerQuantity("1 m ", true).to_string(), string("1 MB"));
    ASSERT_EQ(ComputerQuantity("1 G ", true).to_string(), string("1 GB"));
    ASSERT_EQ(ComputerQuantity("1242Kb", true).to_string(), string("1.21 MB"));
}

TEST(ComputerQuantity, math) {
    ComputerQuantity q{0};
    q += 10;
    ASSERT_EQ( q, 10 );
    q *= 1024;
    ASSERT_EQ( q, 10 * 1024);
    q -= 7;
    ASSERT_EQ( q, 10 * 1024 - 7);

    ComputerQuantity q2{2};
    ComputerQuantity q3{3};
    ComputerQuantity q4 = q2 * q3;
    ASSERT_EQ(q4, 6);

    ASSERT_THROW(operator-(q4, 2048), QuantityError);
}

/**
 * Unit testing for the static method DurationQuantity::parse()
 */
TEST(DurationQuantity, parse){

    // straighforward numbers
    ASSERT_EQ(DurationQuantity::parse("0"), 0 ); // seconds
    ASSERT_EQ(DurationQuantity::parse("7"), 7 ); // 7 seconds
    ASSERT_THROW(DurationQuantity::parse("-1").seconds(), QuantityError );

    // hours
    ASSERT_EQ(DurationQuantity::parse("0h"), 0 );
    ASSERT_EQ(DurationQuantity::parse("1h"),  3600); // 1 hour
    ASSERT_EQ(DurationQuantity::parse("1 hour"),  3600); // 1 hour
    ASSERT_EQ(DurationQuantity::parse("2 hours"),  2 * 3600); // 2 hours
    ASSERT_EQ(DurationQuantity::parse("2h"),  2* 3600); // 2 hours
    ASSERT_EQ(DurationQuantity::parse("2hours"),  2* 3600); // 2 hours
    ASSERT_EQ(DurationQuantity::parse("3 H"),  3* 3600); // 3 hours
    ASSERT_EQ(DurationQuantity::parse("3 Hours"),  3* 3600); // 3 hours
    ASSERT_EQ(DurationQuantity::parse("3 HOUR"),  3* 3600); // 3 hours
    ASSERT_EQ(DurationQuantity::parse("3HoUr"),  3* 3600); // 3 hours
    ASSERT_THROW(DurationQuantity::parse("2 hrs"), QuantityError );

    // minutes
    ASSERT_EQ(DurationQuantity::parse("0m"), 0 );
    ASSERT_EQ(DurationQuantity::parse("1m"),  60); // 1 minute
    ASSERT_EQ(DurationQuantity::parse("1 min"),  60); // 1 minute
    ASSERT_EQ(DurationQuantity::parse("1 minute"),  60); // 1 minute
    ASSERT_EQ(DurationQuantity::parse("2 m"),  2 * 60); // 2 minutes
    ASSERT_EQ(DurationQuantity::parse("2 m"),  2 * 60); // 2 minutes
    ASSERT_EQ(DurationQuantity::parse("2 min"),  2 * 60); // 2 minutes
    ASSERT_EQ(DurationQuantity::parse("2mins"),  2 * 60); // 2 minutes
    ASSERT_EQ(DurationQuantity::parse("2 minutes"),  2 * 60); // 2 minutes
    ASSERT_EQ(DurationQuantity::parse("2  minutes   \t "),  2 * 60); // 2 minutes

    // seconds
    ASSERT_EQ(DurationQuantity::parse("0s"), 0 );
    ASSERT_EQ(DurationQuantity::parse("1 sec"), 1);
    ASSERT_EQ(DurationQuantity::parse("1 second"), 1);
    ASSERT_EQ(DurationQuantity::parse("2 secs"), 2);
    ASSERT_EQ(DurationQuantity::parse("2 seconds"), 2);
    ASSERT_EQ(DurationQuantity::parse("3 s"), 3);

    // milliseconds
    ASSERT_EQ(DurationQuantity::parse("0 ms").as<chrono::milliseconds>(), 0ms ); // milliseconds
    ASSERT_EQ(DurationQuantity::parse("1ms").as<chrono::milliseconds>(), 1ms );
    ASSERT_EQ(DurationQuantity::parse("1 ms").as<chrono::milliseconds>(), 1ms );
    ASSERT_EQ(DurationQuantity::parse("2 millisec").as<chrono::milliseconds>(), 2ms );
    ASSERT_EQ(DurationQuantity::parse("3 millisecs").as<chrono::milliseconds>(), 3ms );
    ASSERT_EQ(DurationQuantity::parse("4 msec").as<chrono::milliseconds>(), 4ms );
    ASSERT_EQ(DurationQuantity::parse("5 msecs").as<chrono::milliseconds>(), 5ms );
    ASSERT_EQ(DurationQuantity::parse("6 millisec").as<chrono::milliseconds>(), 6ms );
    ASSERT_EQ(DurationQuantity::parse("7 milliseconds").as<chrono::milliseconds>(), 7ms );

    // microsecs
    ASSERT_EQ(DurationQuantity::parse("0 us").as<chrono::microseconds>(), 0us ); // microsecs
    ASSERT_EQ(DurationQuantity::parse("1us").as<chrono::microseconds>(), 1us );
    ASSERT_EQ(DurationQuantity::parse("1 us").as<chrono::microseconds>(), 1us );
    ASSERT_EQ(DurationQuantity::parse("2 microsec").as<chrono::microseconds>(), 2us );
    ASSERT_EQ(DurationQuantity::parse("3 microsecs").as<chrono::microseconds>(), 3us );
    ASSERT_EQ(DurationQuantity::parse("4 usec").as<chrono::microseconds>(), 4us );
    ASSERT_EQ(DurationQuantity::parse("5 usecs").as<chrono::microseconds>(), 5us );
    ASSERT_EQ(DurationQuantity::parse("6 microsecond").as<chrono::microseconds>(), 6us );
    ASSERT_EQ(DurationQuantity::parse("7 microseconds").as<chrono::microseconds>(), 7us );

    // nanoseconds
    ASSERT_EQ(DurationQuantity::parse("0 ns").as<chrono::nanoseconds>(), 0ns ); // nanoseconds
    ASSERT_EQ(DurationQuantity::parse("1ns").as<chrono::nanoseconds>(), 1ns );
    ASSERT_EQ(DurationQuantity::parse("1 ns").as<chrono::nanoseconds>(), 1ns );
    ASSERT_EQ(DurationQuantity::parse("2 nanosec").as<chrono::nanoseconds>(), 2ns );
    ASSERT_EQ(DurationQuantity::parse("3 nanosecs").as<chrono::nanoseconds>(), 3ns );
    ASSERT_EQ(DurationQuantity::parse("4 nsec").as<chrono::nanoseconds>(), 4ns );
    ASSERT_EQ(DurationQuantity::parse("5 nsecs").as<chrono::nanoseconds>(), 5ns );
    ASSERT_EQ(DurationQuantity::parse("6 nanosecond").as<chrono::nanoseconds>(), 6ns );
    ASSERT_EQ(DurationQuantity::parse("7 nanoseconds").as<chrono::nanoseconds>(), 7ns );
}
