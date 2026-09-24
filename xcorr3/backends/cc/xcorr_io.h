#ifndef XCORR_IO_H
#define XCORR_IO_H
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>
#include <utility>
#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>

inline std::runtime_error io_error(const std::string &action) {
    return std::runtime_error(action + ": " + std::strerror(errno));
}

class CheckedFile {
    FILE *file_ = nullptr;
    std::string path_;
public:
    explicit CheckedFile(const std::string &path): path_(path) {
        file_ = std::fopen(path.c_str(), "w");
        if (!file_) throw io_error("open " + path);
    }
    CheckedFile(const CheckedFile &) = delete;
    CheckedFile &operator=(const CheckedFile &) = delete;
    ~CheckedFile() { if (file_) std::fclose(file_); }
    FILE *get() const { return file_; }
    void close() {
        if (!file_) return;
        if (std::fflush(file_) != 0 || std::ferror(file_)) throw io_error("flush " + path_);
        int result;
        do { result = fsync(fileno(file_)); } while (result < 0 && errno == EINTR);
        if (result < 0) throw io_error("sync " + path_);
        FILE *closing = file_; file_ = nullptr;
        if (std::fclose(closing) != 0) throw io_error("close " + path_);
    }
};

/* Private same-filesystem staging. No result is published until all requested
 * stages succeed. Each rename is atomic; the group is not a filesystem-wide
 * atomic transaction. Detected publication failures roll back earlier names. */
class OutputTransaction {
    std::string base_, directory_;
    std::vector<std::string> names_;
    bool preserve_ = false;
    static void remove_private(const std::string &path) noexcept {
        DIR *dir = opendir(path.c_str());
        if (!dir) return;
        while (dirent *entry = readdir(dir)) {
            std::string name = entry->d_name;
            if (name == "." || name == "..") continue;
            std::string child = path + "/" + name;
            struct stat st;
            if (lstat(child.c_str(), &st) == 0 && S_ISDIR(st.st_mode)) remove_private(child);
            else unlink(child.c_str()); // Never follow a symlink out of our workspace.
        }
        closedir(dir); rmdir(path.c_str());
    }
    static bool exists_regular(const std::string &path) {
        struct stat st;
        if (lstat(path.c_str(), &st) == 0) {
            if (!S_ISREG(st.st_mode)) throw std::runtime_error("output target is not a regular file: " + path);
            return true;
        }
        if (errno != ENOENT) throw io_error("inspect output " + path);
        return false;
    }
public:
    explicit OutputTransaction(std::vector<std::string> names): names_(std::move(names)) {
        char *cwd = realpath(".", nullptr);
        if (!cwd) throw io_error("resolve working directory");
        base_ = cwd; std::free(cwd);
        for (const auto &name : names_) exists_regular(base_ + "/" + name);
        std::string pattern = base_ + "/.xcorr_cc-XXXXXX";
        std::vector<char> buf(pattern.begin(), pattern.end()); buf.push_back('\0');
        if (!mkdtemp(buf.data())) throw io_error("create output staging directory");
        directory_ = buf.data();
    }
    OutputTransaction(const OutputTransaction &) = delete;
    OutputTransaction &operator=(const OutputTransaction &) = delete;
    ~OutputTransaction() { if (!preserve_) remove_private(directory_); }
    const std::string &directory() const { return directory_; }
    std::string path(const std::string &name) const { return directory_ + "/" + name; }
    static void require_product(const std::string &path) {
        struct stat st;
        if (lstat(path.c_str(), &st) != 0) throw io_error("missing output " + path);
        if (!S_ISREG(st.st_mode) || st.st_size == 0) throw std::runtime_error("empty or nonregular output: " + path);
    }
    void publish() {
        for (const auto &name : names_) require_product(path(name));
        std::vector<bool> backed(names_.size(), false), installed(names_.size(), false);
        try {
            for (size_t i=0; i<names_.size(); ++i) {
                std::string dest = base_ + "/" + names_[i], backup = path("backup-" + std::to_string(i));
                if (exists_regular(dest)) {
                    if (rename(dest.c_str(), backup.c_str()) != 0) throw io_error("backup " + dest);
                    backed[i] = true;
                }
                if (rename(path(names_[i]).c_str(), dest.c_str()) != 0) throw io_error("publish " + dest);
                installed[i] = true;
            }
        } catch (...) {
            bool restored = true;
            for (size_t j=names_.size(); j>0; --j) {
                size_t i=j-1; std::string dest=base_+"/"+names_[i];
                if (installed[i] && unlink(dest.c_str()) != 0) restored=false;
                if (backed[i] && rename(path("backup-"+std::to_string(i)).c_str(),dest.c_str()) != 0) restored=false;
            }
            if (!restored) {
                preserve_ = true;
                throw std::runtime_error("output rollback incomplete; recover saved files from " + directory_);
            }
            throw;
        }
    }
};
#endif
