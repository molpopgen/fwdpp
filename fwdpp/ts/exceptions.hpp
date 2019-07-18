#ifndef FWDPP_TS_EXCEPTIONS_HPP
#define FWDPP_TS_EXCEPTIONS_HPP

#include <string>
#include <stdexcept>

namespace fwdpp
{
    namespace ts
    {
        class tables_error : public std::exception
        /// Thrown when table_collection objects are in an invalid state.
        /// \version 0.8.0 Added to library
        {
          private:
            std::string message_;

          public:
            explicit tables_error(std::string message)
                : message_(std::move(message))
            {
            }
            virtual const char*
            what() const noexcept
            {
                return message_.c_str();
            }
        };

        class samples_error : public std::exception
        /// Thrown when there is a problem with samples lists
        /// \version 0.8.0 Added to library
        {
          private:
            std::string message_;

          public:
            explicit samples_error(std::string message)
                : message_(std::move(message))
            {
            }
            virtual const char*
            what() const noexcept
            {
                return message_.c_str();
            }
        };

    } // namespace ts
} // namespace fwdpp

#endif

