#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

namespace kaspectra::math {

    enum class InterpMode { Linear, LogLog };

    // Edge case: queries outside [xs.font(), xs.back()] is clamped to the
    // boundary y-value rahter than extrapolation. This should be more stable,
    // but can be modified later
    class Interpolator1D {
        public:
            Interpolator1D(std::vector<double> xs, std::vector<double> ys, InterpMode mode = InterpMode::LogLog):
                xs_(std::move(xs)), ys_(std::move(ys)), mode_(mode) {
                if (xs_.size() != ys_.size() || xs_.size() < 2) {
                        throw std::invalid_argument("Interpolator1D: need matching, >=2 point arrays");
                    }
                if (mode_ == InterpMode::LogLog) {
                    for (std::size_t i = 0; i < xs_.size(); ++i) {
                        if (!(xs_[i] > 0.0) || !(ys_[i] > 0.0)) {
                            throw std::invalid_argument("Interpolator1D: LogLog mode requires x>0 and y>0 for all points");
                        }
                    }
                }
            }

            double operator()(double x) const { return at(x); }

            double at(double x) const {
                if (x <= xs_.front()) return ys_.front();
                if (x >= xs_.back()) return ys_.back();

                auto it = std::lower_bound(xs_.begin(), xs_.end(), x);
                std::size_t i = static_cast<std::size_t>(it - xs_.begin());
                std::size_t lo = i - 1, hi = i;

                if (mode_ == InterpMode::Linear) {
                    double t = (x - xs_[lo]) / (xs_[hi] - xs_[lo]);
                    return ys_[lo] + t * (ys_[hi] - ys_[lo]);
                }

                double lx = std::log(x);
                double lx0 = std::log(xs_[lo]), lx1 = std::log(xs_[hi]);
                double ly0 = std::log(ys_[lo]), ly1 = std::log(ys_[hi]);
                double t = (lx - lx0) / (lx1 - lx0);
                return std::exp(ly0 + t * (ly1 - ly0));
            }
        
        private:
            std::vector<double> xs_, ys_;
            InterpMode mode_;
    };

// Usage:
//   Interpolator1D phi({1e-6, 1e-3, 1.0}, {2.3e-4, 1.1e-2, 0.0}, InterpMode::LogLog);
//   double y = phi(5e-4)

}  // namespace kaspectra::math