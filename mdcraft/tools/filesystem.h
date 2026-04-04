#pragma once

#if __has_include(<filesystem>)
#include <filesystem>
namespace mdcraft { namespace tools { namespace filesystem = std::filesystem; } }
#elif __has_include(<experimental/filesystem>)
#include <experimental/filesystem>
namespace mdcraft { namespace tools { namespace filesystem = std::experimental::filesystem; } }
#else
#error "Neither <filesystem> nor <experimental/filesystem> is available"
#endif
