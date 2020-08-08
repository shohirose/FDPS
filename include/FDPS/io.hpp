#pragma once

#include <string>

namespace ParticleSimulator {

/// @brief Checks if a given file or directory exist.
/// @param[in] fileOrDirectoryName A file or directory name
/// @returns True if the given file or directory exists.
bool exists(const std::string& fileOrDirectoryName) noexcept;

/// @brief Creates a directory of a given name.
/// @param[in] directoryName A directory name
/// @returns True if a directory is created successfully.
bool createDirectory(const std::string& directoryName) noexcept;

/// @brief Change permissions of a given file or directory.
/// @param[in] fileOrDirectoryName A file or directory name
/// @param[in] mode Permission mode in the format of UNIX's 4-digits octal
/// number.
/// @returns True if the permission is changed successfully.
bool permissions(const std::string& fileOrDirectoryName, int mode) noexcept;

}  // namespace ParticleSimulator
