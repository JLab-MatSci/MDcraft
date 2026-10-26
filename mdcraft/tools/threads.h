#pragma once

#include <vector>
#include <limits>
#include <algorithm>

#include <mdcraft/configuration.h>

#ifdef mdcraft_ENABLE_TBB
#include <execution>
#include <tbb/global_control.h>
#include <tbb/parallel_for.h>
#else
#include <mdcraft/tools/thread-pool.h>
#endif

#ifdef mdcraft_ENABLE_MPI
#include <mdcraft/tools/network.h>
#endif

namespace mdcraft::tools {

/// @brief Number of tasks per thread by default (using thread-pool)
constexpr int default_n_tasks_per_thread = 10;


/// @brief Simplified interface for multi-threading
class Threads {
public:
    /// @param count Number of threads. Large count will be limited by the
    /// number of cores. Use count < 2 to disable multi-threading.
    explicit Threads(int count = 1000000);

#ifdef mdcraft_ENABLE_MPI
    /// @param count Number of threads. Large count will be limited by
    /// the number of cores. Use count < 2 to disable multi-threading.
    /// @param comm MPI-communicator is used to determine an optimal
    /// number of threads.
    Threads(int count, MPI_Comm comm);
#endif

    /// @brief Enable multi-threading using all available cores
    void on();

    /// @brief Disable multi-threading
    void off();

    /// @brief Number of threads
    int count() { return m_count; }

    /// @brief Maximum number of threads, all available cores,
    /// taking into account the number of MPI processes
    int max_count() { return m_max_count; }

    /// @brief Multi-threading is active?
    bool active() { return m_count > 1; }

    /// @brief Multi-threading is disabled?
    bool disabled() { return m_count < 2; }

#ifndef mdcraft_ENABLE_TBB
    ThreadPool& pool() { return *m_pool; }
#endif

    /// @brief Execute function for elements in range [begin, end)
    /// @param begin, end Iterators
    /// @param func Target function, take arguments (*Iter, Args...)
    /// @param args Arguments of target function
    /// @details The target function takes as arguments a dereferenced
    /// iterator *Iter and a pack of arguments Args..., the return value
    /// of target function is ignored.
    /// Короче функция может использовать просто итераторы func(Iter, Args...),
    /// я добавил такую возможность, потому что половина кода на итераторах,
    /// но вообще говоря это не круто!!!
    template<int n_tasks_per_thread = default_n_tasks_per_thread,
            typename Iter, class Func, class... Args>
    void for_each(Iter begin, Iter end, Func &&func, Args &&... args);

    /// @brief Execute function for indices in range [begin, end)
    /// @param func Target function, take arguments (Index, Args...)
    /// @param args Arguments of target function
    /// @details The target function takes as arguments an index and a pack
    /// of arguments Args..., the return value of target function is ignored.
    template<int n_tasks_per_thread = default_n_tasks_per_thread,
            typename Index, class Func, class... Args,
            typename = std::enable_if_t<std::is_integral_v<Index>>>
    void parallel_for(Index begin, Index end, Func &&func, Args &&... args);

    /// @brief Minimize the results of executing a function on a range of
    /// elements [begin, end).
    /// @param begin, end Iterators
    /// @param init Initial "largest" value
    /// @param func Target function, take arguments (*Iter, Args...)
    /// @param args Arguments of target function
    /// @tparam Value The return type of the target function; the comparison
    /// operator< "less" must be defined for this type.
    /// @details The target function takes as arguments a dereferenced
    /// iterator *Iter and a pack of arguments Args..., the return value
    /// of target function has type Value.
    template<int n_tasks_per_thread = default_n_tasks_per_thread,
            typename Iter, typename Value, class Func, class... Args>
    Value min(Iter begin, Iter end, const Value &init, Func &&func, Args &&... args);

    /// @brief Minimize the results of executing a function on a range of
    /// elements [begin, end). The return type of the target function must
    /// be an arithmetic type.
    /// @param begin, end Iterators
    /// @param func Target function, take arguments (*Iter, Args...)
    /// @param args Arguments of target function
    /// @details The target function takes as arguments a dereferenced
    /// iterator *Iter and a pack of arguments Args...
    template<int n_tasks_per_thread = default_n_tasks_per_thread,
            typename Iter, class Func, class... Args>
    auto min(Iter begin, Iter end, Func &&func, Args &&... args)
    -> std::invoke_result_t<Func, decltype(*begin), Args...> {
        using Value = std::invoke_result_t<Func, decltype(*begin), Args...>;
        static_assert(std::is_arithmetic_v<Value>, "Function must have arithmetic return type");
        return min<n_tasks_per_thread>(begin, end, std::numeric_limits<Value>::max(),
                                       std::forward<Func>(func), std::forward<Args>(args)...);
    }

    /// @brief Maximize the results of executing a function on a range of
    /// elements [begin, end).
    /// @param begin, end Iterators
    /// @param init Initial "lowest" value
    /// @param func Target function, take arguments (*Iter, Args...)
    /// @param args Arguments of target function
    /// @tparam Value The return type of the target function; the comparison
    /// operator> "greater" must be defined for this type.
    /// @details The target function takes as arguments a dereferenced
    /// iterator *Iter and a pack of arguments Args..., the return value
    /// of target function has type Value.
    template<int n_tasks_per_thread = default_n_tasks_per_thread,
            typename Iter, typename Value, class Func, class... Args>
    Value max(Iter begin, Iter end, const Value &init, Func &&func, Args &&... args);

    /// @brief Maximize the results of executing a function on a range of
    /// elements [begin, end). The return type of the target function must
    /// be an arithmetic type.
    /// @param begin, end Iterators
    /// @param func Target function, take arguments (*Iter, Args...)
    /// @param args Arguments of target function
    /// @details The target function takes as arguments a dereferenced
    /// iterator *Iter and a pack of arguments Args...
    template<int n_tasks_per_thread = default_n_tasks_per_thread,
            typename Iter, class Func, class... Args>
    auto max(Iter begin, Iter end, Func &&func, Args &&... args)
    -> std::invoke_result_t<Func, decltype(*begin), Args...> {
        using Value = std::invoke_result_t<Func, decltype(*begin), Args...>;
        static_assert(std::is_arithmetic_v<Value>, "Function must have arithmetic return type");
        return max<n_tasks_per_thread>(begin, end, std::numeric_limits<Value>::lowest(),
                                       std::forward<Func>(func), std::forward<Args>(args)...);
    }

    /// @brief Summarize the results of executing a function on a range of
    /// elements [begin, end).
    /// @param begin, end Iterators
    /// @param init Initial (zero) value
    /// @param func Target function, take arguments (*Iter, Args...)
    /// @param args Arguments of target function
    /// @tparam Value The return type of the target function; the append
    /// operator+= must be defined for this type.
    /// @details The target function takes as arguments a dereferenced
    /// iterator *Iter and a pack of arguments Args..., the return value
    /// of target function has type Value.
    template<int n_tasks_per_thread = default_n_tasks_per_thread,
            typename Iter, typename Value, class Func, class... Args>
    Value sum(Iter begin, Iter end, const Value &init, Func &&func, Args &&... args);

    /// @brief Common reducing operation on the results of executing a function
    /// on a range of elements [begin, end). Using operator &=.
    /// @param begin, end Iterators
    /// @param init Initial value
    /// @param func Target function, take arguments (*Iter, Args...)
    /// @param args Arguments of target function
    /// @tparam Value The return type of the target function; the operator&=
    /// must be defined for this type.
    /// @details The target function takes as arguments a dereferenced
    /// iterator *Iter and a pack of arguments Args..., the return value
    /// of target function has type Value.
    template<int n_tasks_per_thread = default_n_tasks_per_thread,
            typename Iter, typename Value, class Func, class... Args>
    Value reduce(Iter begin, Iter end, const Value &init, Func &&func, Args &&... args);

    /// @brief Parallel sort
    /// @param begin, end Iterators
    /// @tparam BinaryOp Comparator
    template<typename Iter, typename BinaryOp = std::less<typename Iter::value_type>>
    void sort(const Iter& begin, const Iter& end, BinaryOp comp = BinaryOp());

private:
    /// @brief Threads count (< 2 -- serial execution)
    int m_count = 1;

    /// @brief Recommended number of threads, all available cores,
    /// taking into account the number of MPI processes
    int m_max_count;

#ifdef mdcraft_ENABLE_TBB
    /// @brief Using TBB parallelism
    std::shared_ptr<tbb::global_control> m_control = nullptr;
#else
    /// @brief Using simple thread pool
    std::shared_ptr<ThreadPool> m_pool = nullptr;
#endif


    /// @brief Enable multi-threading
    /// @param count Number of threads
    void resize(int count);

    /// @brief Optimal number of task for the thread pool
    /// @param n_elements Number of elements in range
    template<int n_tasks_per_thread = default_n_tasks_per_thread>
    int get_n_tasks(size_t n_elements) const {
        // Minimal number of elements per task (using thread-pool)
        constexpr int min_elements_per_task = 500;
        return std::max(1, std::min(n_tasks_per_thread * m_count,
                                    int(n_elements / min_elements_per_task)));
    }
};

extern Threads dummy_pool;

template<int n_tpt, typename Iter, class Func, class ...Args>
void Threads::for_each(Iter begin, Iter end, Func &&func, Args &&... args) {
    constexpr bool func_use_iter = std::is_invocable_v<std::decay_t<Func>, Iter, Args...>;

    if (Threads::disabled()) { // serial execution
        if constexpr (func_use_iter) {
            // if function use iterator
            for (auto it = begin; it < end; ++it) {
                func( it, std::forward<Args>(args)...);
            }
        }
        else {
            // if function use dereferenced value
            for (auto it = begin; it < end; ++it) {
                func(*it, std::forward<Args>(args)...);
            }
        }
        return;
    }

#ifdef mdcraft_ENABLE_TBB
    if constexpr (func_use_iter) {
        // if function use iterator
        tbb::parallel_for(
                      size_t{0}, static_cast<size_t>(end - begin),
                      [begin, &func, &args...](size_t idx) {
                          func(begin + idx, std::forward<Args>(args)...);
                      });
    }
    else {
        // if function use dereferenced value
        std::for_each(std::execution::par,
                      begin, end,
                      [&func, &args...](auto &&elem) {
                          func(elem, std::forward<Args>(args)...);
                      });
    }
    return;
#else
    auto bin_function =
            [&func, &args...](const Iter &a, const Iter &b) {
                if constexpr (func_use_iter) {
                    // if function use iterator
                    for (auto it = a; it < b; ++it) {
                        func( it, std::forward<Args>(args)...);
                    }
                }
                else {
                    // if function use dereferenced value
                    for (auto it = a; it < b; ++it) {
                        func(*it, std::forward<Args>(args)...);
                    }
                }
            };

    size_t size = end - begin;

    // empty range
    if (size < 1) return;

    int n_tasks = get_n_tasks<n_tpt>(size);
    size_t bin = size / n_tasks;
    std::vector<std::future<void>> results;
    results.reserve(n_tasks);

    Iter from = begin;
    for (int i = 0; i < n_tasks - 1; ++i) {
        results.emplace_back(m_pool->enqueue(bin_function, from, from + bin));
        from += bin;
    }
    results.emplace_back(m_pool->enqueue(bin_function, from, end));

    for (auto &result : results) result.wait();
    for (auto &result: results)
        result.get();
    return;
#endif
}

template<int n_tpt, typename Index, class Func, class ...Args, typename >
void Threads::parallel_for(Index begin, Index end, Func &&func, Args &&... args) {
    if (Threads::disabled()) { // serial execution
        for (Index i = begin; i < end; ++i) {
            func(i, std::forward<Args>(args)...);
        }
        return;
    }

#ifdef mdcraft_ENABLE_TBB
    tbb::parallel_for(begin, end,
                      [&func, &args...](Index idx) {
                          func(idx, std::forward<Args>(args)...);
                      });
    return;
#else
    auto bin_function =
            [&func, &args...](Index from, Index to) {
                for (Index idx = from; idx < to; ++idx) {
                    func(idx, std::forward<Args>(args)...);
                }
            };

    Index size = end - begin;

    // empty range
    if (size < 1) return;

    int n_tasks = get_n_tasks<n_tpt>(size);
    size_t bin = size / n_tasks;
    std::vector<std::future<void>> results;
    results.reserve(n_tasks);

    Index from = begin;
    for (int i = 0; i < n_tasks - 1; ++i) {
        results.emplace_back(m_pool->enqueue(bin_function, from, from + bin));
        from += bin;
    }
    results.emplace_back(m_pool->enqueue(bin_function, from, end));

    for (auto &result : results) result.wait();
    for (auto &result : results)
        result.get();
    return;
#endif
}

template<int n_tpt, typename Iter, typename Value, class Func, class ...Args>
Value Threads::min(Iter begin, Iter end, const Value &init, Func &&func, Args &&... args) {
    using ReturnType = std::invoke_result_t<Func, decltype(*begin), Args...>;
    static_assert(std::is_convertible_v<ReturnType, Value>,
                  "Function return type and initial value type is not consistent");

    if (Threads::disabled()) { // serial execution
        Value res{init};
        for (auto it = begin; it < end; ++it) {
            Value temp = func(*it, std::forward<Args>(args)...);
            if (temp < res) { res = temp; }
        }
        return res;
    }

#ifdef mdcraft_ENABLE_TBB
    return std::transform_reduce(std::execution::par,
                                 begin, end, init,
                                 [](auto&& a, auto&& b) { return a < b ? a : b; },
                                 [&func, &args...](auto&& elem) -> Value {
                                     return func(elem, std::forward<Args>(args)...);
                                 });
#else
    auto bin_function =
            [&init, &func, &args...](const Iter &a, const Iter &b) -> Value {
                Value res{init};
                Value temp{init};
                for (auto it = a; it < b; ++it) {
                    temp = func(*it, std::forward<Args>(args)...);
                    if (temp < res) { res = temp; }
                }
                return res;
            };

    int size = end - begin;

    // empty range
    if (size < 1) return init;

    int n_tasks = get_n_tasks<n_tpt>(size);
    size_t bin = size / n_tasks;
    std::vector<std::future<Value>> results;
    results.reserve(n_tasks);

    Iter from = begin;
    for (int i = 0; i < n_tasks - 1; ++i) {
        results.emplace_back(m_pool->enqueue(bin_function, from, from + bin));
        from += bin;
    }
    results.emplace_back(m_pool->enqueue(bin_function, from, end));

    Value res{init};
    Value temp{init};
    for (auto &result : results) {
        temp = result.get();
        if (temp < res) { res = temp; }
    }
    return res;
#endif
}

template<int n_tpt, typename Iter, typename Value, class Func, class ...Args>
Value Threads::max(Iter begin, Iter end, const Value &init, Func &&func, Args &&... args) {
    using ReturnType = std::invoke_result_t<Func, decltype(*begin), Args...>;
    static_assert(std::is_convertible_v<ReturnType, Value>,
                  "Function return type and initial value type is not consistent");

    if (Threads::disabled()) { // serial execution
        Value res{init};
        for (auto it = begin; it < end; ++it) {
            Value temp = func(*it, std::forward<Args>(args)...);
            if (temp > res) { res = temp; }
        }
        return res;
    }

#ifdef mdcraft_ENABLE_TBB
    return std::transform_reduce(std::execution::par,
                                 begin, end, init,
                                 [](auto&& a, auto&& b) { return a > b ? a : b; },
                                 [&func, &args...](auto&& elem) -> Value {
                                     return func(elem, std::forward<Args>(args)...);
                                 });
#else
    auto bin_function =
            [&init, &func, &args...](const Iter &a, const Iter &b) -> Value {
                Value res{init};
                Value temp{init};
                for (auto it = a; it < b; ++it) {
                    temp = func(*it, std::forward<Args>(args)...);
                    if (temp > res) { res = temp; }
                }
                return res;
            };

    int size = end - begin;

    // empty range
    if (size < 1) return init;

    int n_tasks = get_n_tasks<n_tpt>(size);
    size_t bin = size / n_tasks;
    std::vector<std::future<Value>> results;
    results.reserve(n_tasks);

    Iter from = begin;
    for (int i = 0; i < n_tasks - 1; ++i) {
        results.emplace_back(m_pool->enqueue(bin_function, from, from + bin));
        from += bin;
    }
    results.emplace_back(m_pool->enqueue(bin_function, from, end));

    Value res{init};
    Value temp{init};
    for (auto &result : results) {
        temp = result.get();
        if (temp > res) { res = temp; }
    }
    return res;
#endif
}

template<int n_tpt, typename Iter, typename Value, class Func, class ...Args>
Value Threads::sum(Iter begin, Iter end, const Value &init, Func &&func, Args &&... args) {
    using ReturnType = std::invoke_result_t<Func, decltype(*begin), Args...>;
    static_assert(std::is_convertible_v<ReturnType, Value>,
                  "Function return type and initial value type is not consistent");

    if (Threads::disabled()) { // serial execution
        Value res{init};
        for (auto it = begin; it < end; ++it) {
            res += func(*it, std::forward<Args>(args)...);
        }
        return res;
    }

#ifdef mdcraft_ENABLE_TBB
    return std::transform_reduce(std::execution::par,
                                 begin, end, init,
                                 [](auto&& a, auto&& b) { auto c = a; c += b; return c; },
                                 [&func, &args...](auto&& elem) -> Value {
                                     return func(elem, std::forward<Args>(args)...);
                                 });
#else
    auto bin_function =
            [&init, &func, &args...](const Iter &a, const Iter &b) -> Value {
                Value res{init};
                for (auto it = a; it < b; ++it) {
                    res += func(*it, std::forward<Args>(args)...);
                }
                return res;
            };

    int size = end - begin;

    // empty range
    if (size < 1) return init;

    int n_tasks = get_n_tasks<n_tpt>(size);
    size_t bin = size / n_tasks;
    std::vector<std::future<Value>> results;
    results.reserve(n_tasks);

    Iter from = begin;
    for (int i = 0; i < n_tasks - 1; ++i) {
        results.emplace_back(m_pool->enqueue(bin_function, from, from + bin));
        from += bin;
    }
    results.emplace_back(m_pool->enqueue(bin_function, from, end));

    Value res{init};
    for (auto &result : results) {
        res += result.get();
    }
    return res;
#endif
}

template<int n_tpt, typename Iter, typename Value, class Func, class ...Args>
Value Threads::reduce(Iter begin, Iter end, const Value &init, Func &&func, Args &&... args) {
    using ReturnType = std::invoke_result_t<Func, decltype(*begin), Args...>;
    static_assert(std::is_convertible_v<ReturnType, Value>,
                  "Function return type and initial value type is not consistent");

    if (Threads::disabled()) { // serial execution
        Value res{init};
        for (auto it = begin; it < end; ++it) {
            res &= func(*it, std::forward<Args>(args)...);
        }
        return res;
    }

#ifdef mdcraft_ENABLE_TBB
    return std::transform_reduce(std::execution::par,
                                 begin, end, init,
                                 [](auto&& a, auto&& b) { auto c = a; c &= b; return c; },
                                 [&func, &args...](auto&& elem) -> Value {
                                     return func(elem, std::forward<Args>(args)...);
                                 });
#else
    auto bin_function =
            [&init, &func, &args...](const Iter &a, const Iter &b) -> Value {
                Value res{init};
                for (auto it = a; it < b; ++it) {
                    res &= func(*it, std::forward<Args>(args)...);
                }
                return res;
            };

    int size = end - begin;

    // empty range
    if (size < 1) return init;

    int n_tasks = get_n_tasks<n_tpt>(size);
    size_t bin = size / n_tasks;
    std::vector<std::future<Value>> results;
    results.reserve(n_tasks);

    Iter from = begin;
    for (int i = 0; i < n_tasks - 1; ++i) {
        results.emplace_back(m_pool->enqueue(bin_function, from, from + bin));
        from += bin;
    }
    results.emplace_back(m_pool->enqueue(bin_function, from, end));

    Value res{init};
    for (auto &result : results) {
        res &= result.get();
    }
    return res;
#endif
}

template<typename Iter, typename BinaryOp = std::less<typename Iter::value_type>>
std::vector<typename Iter::value_type>
merge(const Iter& l, const Iter& m, const Iter& r, BinaryOp u = BinaryOp()) {
	auto size = std::distance(l, r);
	std::vector<typename Iter::value_type> result;
	result.reserve(size);

	auto i = l;
	auto j = m;

	while (i != m && j != r) {
		auto b = u(*i, *j);
		result.push_back(b || *i == *j ? *i++ : *j++);
	}

	result.insert(result.end(), i, m);
	result.insert(result.end(), j, r);

	return result;
}

template<typename Iter, typename BinaryOp>
void Threads::sort(const Iter& begin, const Iter& end, BinaryOp comp) {
    if (Threads::disabled()) { // serial execution
        std::sort(begin, end, comp);
        return;
    }

#ifdef mdcraft_ENABLE_TBB
    std::sort(std::execution::par, begin, end, comp);
    return;
#else
    // Merge sorting using thread pool
	std::size_t size = std::distance(begin, end);
	std::size_t num_tasks = m_pool->size();
	std::size_t one_task_size = std::min(size / num_tasks + 1, size);
	num_tasks = std::min(size / one_task_size + 1, num_tasks);

	{
		std::vector<std::future<void>> results;

		auto i1 = begin;
		auto i2 = begin + one_task_size;

		for (int i = 0; i < num_tasks; i++) {

			auto lambda = [](Iter i, Iter j, BinaryOp u){ std::sort(i, j, u); };

			std::future<void> res = m_pool->enqueue(std::move(lambda), i1, i2, comp);
			results.push_back(std::move(res));

			i1 += one_task_size;
			i2 += one_task_size;
			if (i2 > end) i2 = end;
		}
		for (auto &res : results) res.wait();
		std::for_each(results.begin(), results.end(), [](std::future<void>& res) {res.get();});
	}

	{
		using res_type = std::vector<typename Iter::value_type>;
		while (one_task_size < size) {
			std::vector<std::future<res_type>> results;
			Iter i = begin;
			while (i < end) {
				auto i1 = i;
				auto i2 = std::min(i + one_task_size, end);
				auto i3 = std::min(i + 2. * one_task_size, end);

				std::future<res_type> res = m_pool->enqueue([i1, i2, i3, comp](){return merge(i1,i2,i3,comp);});
				results.push_back(std::move(res));

				i += 2 * one_task_size;
			}
			i = begin;
			std::for_each(results.begin(), results.end(), [&](std::future<res_type>& res) {
				auto&& vec = res.get();
				std::move(vec.begin(), vec.end(), i);
				i += 2 * one_task_size;
			});

			one_task_size *= 2;
		}
	}
#endif
}

} // namespace mdcraft::tools
