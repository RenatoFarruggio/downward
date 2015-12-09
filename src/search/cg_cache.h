#ifndef CG_CACHE_H
#define CG_CACHE_H

#include "task_proxy.h"

#include <vector>

struct ValueTransitionLabel;

class CGCache {
    TaskProxy task_proxy;
    std::vector<std::vector<int>> cache;
    std::vector<std::vector<ValueTransitionLabel *>> helpful_transition_cache;
    std::vector<std::vector<int>> depends_on;

    int get_index(int var, const State &state, int from_val, int to_val) const;
    int compute_required_cache_size(int var_id,
                                    const std::vector<int> &depends_on) const;
public:
    static const int NOT_COMPUTED = -2;

    explicit CGCache(TaskProxy &task_proxy);
    ~CGCache();

    bool is_cached(int var) const {
        return !cache[var].empty();
    }

    int lookup(int var, const State &state, int from_val, int to_val) const {
        return cache[var][get_index(var, state, from_val, to_val)];
    }

    void store(int var, const State &state,
               int from_val, int to_val, int cost) {
        cache[var][get_index(var, state, from_val, to_val)] = cost;
    }

    ValueTransitionLabel *lookup_helpful_transition(
        int var, const State &state, int from_val, int to_val) const {
        int index = get_index(var, state, from_val, to_val);
        return helpful_transition_cache[var][index];
    }

    void store_helpful_transition(
        int var, const State &state, int from_val, int to_val,
        ValueTransitionLabel *helpful_transition) {
        int index = get_index(var, state, from_val, to_val);
        helpful_transition_cache[var][index] = helpful_transition;
    }
};

#endif
