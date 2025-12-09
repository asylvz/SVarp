#ifndef __LOGFILE_H
#define __LOGFILE_H

#include <string>
#include <iosfwd>
#include <ostream>

class LogFile
{
public:
    LogFile();
    ~LogFile();

    bool open(const std::string& path);
    void close();
    bool is_open() const;
    std::ostream& out();

    // convenience builder so you can write params.fp_logs << "text";
    template<typename T>
    LogFile& operator<<(const T& v) {
        out() << v;
        return *this;
    }
    // manipulator (e.g. std::endl)
    LogFile& operator<<(std::ostream& (*pf)(std::ostream&)) {
        out() << pf;
        return *this;
    }

    // Disable copying
    LogFile(const LogFile&) = delete;
    LogFile& operator=(const LogFile&) = delete;

    // Move allowed
    LogFile(LogFile&&) noexcept;
    LogFile& operator=(LogFile&&) noexcept;

private:
    struct Impl;
    Impl* impl;
};

#endif
