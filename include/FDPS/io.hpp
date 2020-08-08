#pragma once

#include <iostream>
#include <string>

// Use std::filesystem
#ifdef __cpp_lib_filesystem
#include <filesystem>
// Use MSVC functions
#elif defined(_MSC_VER)
#include <direct.h>    // _mkdir
#include <errno.h>     // errno
#include <io.h>        // _access, _chmod
#include <sys/stat.h>  // stat
#include <sys/types.h>

#include <cstdio>
#include <cstdlib>
// For GCC or Clang
#else
#include <errno.h>     // errno
#include <sys/stat.h>  // stat
#include <sys/types.h>
#endif

namespace ParticleSimulator {

bool exists(const std::string& fileOrDirectoryName) noexcept {
#ifdef __cpp_lib_filesystem
  namespace fs = std::filesystem;
  return fs::exists(fs::path(fileOrDirectoryName));
#elif defined(_MSC_VER)
  return _access(fileOrDirectoryName.c_str(), 0) == 0;
#else
  struct stat buf;
  if (stat(fileOrDirectoryName.c_str(), &buf) == 0) {
    return S_ISREG(buf.st_mode) || S_ISDIR(buf.st_mode);
  } else {
    switch (errno) {
      case ENOENT:
        std::cerr << "Does not exist: " << fileOrDirectoryName << std::endl;
        break;
      case EACCESS:
        std::cerr << "Could not access to: " << fileOrDirectoryName << '.'
                  << std::endl;
        break;
      case ENAMETOOLONG:
        std::cerr << "Name is too long: " << fileOrDirectoryName << std::endl;
        break;
      case ENOTDIR:
        std::cerr << "Incorrect path name: " << fileOrDirectoryName
                  << std::endl;
      default:
        std::cerr << "Error occured when checking: " << fileOrDirectoryName
                  << "\nError code: " << errno << std::endl;
    }
    return false;
  }
#endif
}

bool createDirectory(const std::string& directoryName) noexcept {
#ifdef __cpp_lib_filesystem
  try {
    namespace fs = std::filesystem;
    return fs::create_directory(fs::path(directoryName));
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
    return false;
  }
#elif _MSC_VER
  return _mkdir(directoryName.c_str()) == 0;
#else
  return mkdir(directoryName.c_str(), 0777);
#endif
}

bool permissions(const std::string& fileOrDirectoryName, int mode) noexcept {
#ifdef __cpp_lib_filesystem
  try {
    namespace fs = std::filesystem;

    fs::perms fsmode = fs::perms::none;
    switch (mode) {
      case 0:
        break;
      case 0400:
        fsmode = fs::perms::owner_read;
        break;
      case 0200:
        fsmode = fs::perms::owner_write;
        break;
      case 0100:
        fsmode = fs::perms::owner_exec;
        break;
      case 0700:
        fsmode = fs::perms::owner_all;
        break;
      case 040:
        fsmode = fs::perms::group_read;
        break;
      case 020:
        fsmode = fs::perms::group_write;
        break;
      case 010:
        fsmode = fs::perms::group_exec;
        break;
      case 070:
        fsmode = fs::perms::group_all;
        break;
      case 04:
        fsmode = fs::perms::others_read;
        break;
      case 02:
        fsmode = fs::perms::others_write;
        break;
      case 01:
        fsmode = fs::perms::others_exec;
        break;
      case 07:
        fsmode = fs::perms::others_all;
        break;
      case 0777:
        fsmode = fs::perms::all;
        break;
      case 04000:
        fsmode = fs::perms::set_uid;
        break;
      case 02000:
        fsmode = fs::perms::set_gid;
        break;
      case 01000:
        fsmode = fs::perms::sticky_bit;
        break;
      case 07777:
        fsmode = fs::perms::mask;
        break;
      default:
        std::cerr << "Unknown permission mode: " << mode << std::endl;
        return false;
    }

    fs::permissions(fs::path(fileOrDirectoryName), fsmode);
    return true;
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
    return false;
  }
#elif defined(_MSC_VER)
  const auto hasSucceeded = _chmod(fileOrDirectoryName.c_str(), _S_IWRITE) == 0;
  if (!hasSucceeded) {
    switch (errno) {
      case ENOENT:
        std::cerr << "Could not find: " << fileOrDirectoryName << std::endl;
        break;
      case EINVAL:
        std::cerr << "Invalid parameter: " << fileOrDirectoryName << std::endl;
        break;
      default:
        std::cerr << "Unknown error occured when changing permission of: "
                  << fileOrDirectoryName << "\nError code: " << errno
                  << std::endl;
    }
  }
  return hasSucceeded;
#else
  const auto hasSucceeded = chmod(fileOrDirectoryName.c_str(), mode) == 0);
  if (!hasSucceeded) {
    switch (errno) {
      case ENOENT:
        std::cerr << "Does not exist: " << fileOrDirectoryName << std::endl;
        break;
      case EACCESS:
        std::cerr << "Could not access to: " << fileOrDirectoryName << '.'
                  << std::endl;
        break;
      case ENAMETOOLONG:
        std::cerr << "Name is too long: " << fileOrDirectoryName << std::endl;
        break;
      case ENOTDIR:
        std::cerr << "Incorrect path name: " << fileOrDirectoryName
                  << std::endl;
      default:
        std::cerr << "Error occured when changing the permission of: "
                  << fileOrDirectoryName << "\nError code: " << errno
                  << std::endl;
    }
  }
  return hasSucceeded;
#endif
}

}  // namespace ParticleSimulator
