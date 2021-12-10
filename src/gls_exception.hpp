//
// Created by Fabio Riccardi on 12/9/21.
//

#ifndef GLS_EXCEPTION_HPP
#define GLS_EXCEPTION_HPP

#include <exception>

namespace gls {

    class exception : public std::exception {
        const std::string _message;

    public:
        explicit exception(std::string message) : _message(std::move(message)) {}

        [[nodiscard]] const char *what() const noexcept override { return _message.c_str(); }
    };

} // namespace gls

#endif  // GLS_EXCEPTION_HPP
