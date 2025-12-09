#include "logfile.h"
#include <fstream>

struct LogFile::Impl {
    std::ofstream ofs;
    bool open(const std::string& path) {
        ofs.open(path);
        return ofs.is_open();
    }
    void close() {
        if (ofs.is_open()) ofs.close();
    }
    bool is_open() const { return ofs.is_open(); }
    std::ostream& out() { return ofs; }
};

LogFile::LogFile(): impl(new Impl()) {}

LogFile::~LogFile() {
    if (impl) {
        impl->close();
        delete impl;
        impl = nullptr;
    }
}

bool LogFile::open(const std::string& path) {
    return impl && impl->open(path);
}
void LogFile::close() { if (impl) impl->close(); }

bool LogFile::is_open() const { return impl && impl->is_open(); }

std::ostream& LogFile::out() { return impl->out(); }

LogFile::LogFile(LogFile&& other) noexcept : impl(other.impl) { other.impl = nullptr; }
LogFile& LogFile::operator=(LogFile&& other) noexcept { if (this != &other) { delete impl; impl = other.impl; other.impl = nullptr; } return *this; }
