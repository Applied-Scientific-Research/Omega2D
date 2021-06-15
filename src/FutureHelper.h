/*
 * FutureHelper.h - non-class methods for futures
 *
 * (c)2021 Applied Scientific Research, Inc.
 *         Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include <future>
#include <chrono>

template <class T>
bool is_future_ready(std::future<T> const& f) {
    if (!f.valid()) return false;
    return f.wait_for(std::chrono::seconds(0)) == std::future_status::ready;
    //return f.wait_for(std::chrono::milliseconds(1)) == std::future_status::ready;
}

