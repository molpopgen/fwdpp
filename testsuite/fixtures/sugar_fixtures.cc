#include "sugar_fixtures.hpp"

// MULTILOCUS
const int mlocuspop_popgenmut_fixture::nloci = 4;

std::vector<mlocuspop_popgenmut_fixture::mutmodel>
mlocuspop_popgenmut_fixture::make_mutmodels()
{
    std::vector<mlocuspop_popgenmut_fixture::mutmodel> rv;
    for (int i = 0; i < nloci; ++i)
        {
            auto mm =
                [this,
                 i](std::queue<std::size_t> &recbin,
                    mlocuspop_popgenmut_fixture::poptype::mcont_t &mutations) {
                    return fwdpp::infsites_popgenmut(
                        recbin, mutations, rng.get(), pop.mut_lookup,
                        this->generation, 0.5,
                        [this, i]() {
                            return gsl_ran_flat(rng.get(), i, i + 1);
                        },
                        []() { return -0.01; }, []() { return 1.; });
                };
            rv.push_back(std::move(mm));
        }
    return rv;
}

std::vector<mlocuspop_popgenmut_fixture::recmodel>
mlocuspop_popgenmut_fixture::make_recmodels()
{
    std::vector<mlocuspop_popgenmut_fixture::recmodel> rv;
    for (int i = 0; i < nloci; ++i)
        {
            rv.emplace_back(fwdpp::recbinder(
                fwdpp::poisson_xover(5e-3, i, i + 1), this->rng.get()));
        }
    return rv;
}

std::vector<fwdpp::extensions::discrete_rec_model>
mlocuspop_popgenmut_fixture::fill_vdrm()
{
    double recrate_region = 1e-4;
    // set up a vector of extensions::discrete_rec_model
    // with different regions sizes and weights
    std::vector<fwdpp::extensions::discrete_rec_model> vdrm_;
    for (int i = 0; i < nloci; ++i)
        {
            std::vector<fwdpp::extensions::discrete_rec_model::function_type>
                f;
            std::vector<double> w{ 1., 10., 1. };
            f.push_back([this, i](std::vector<double> &b) {
                b.push_back(gsl_ran_flat(rng.get(), static_cast<double>(i),
                                         static_cast<double>(i) + 1. / 3.));
            });
            f.push_back([this, i](std::vector<double> &b) {
                b.push_back(gsl_ran_flat(rng.get(),
                                         static_cast<double>(i) + 1. / 3,
                                         static_cast<double>(i) + 2. / 3.));
            });
            f.push_back([this, i](std::vector<double> &b) {
                b.push_back(gsl_ran_flat(rng.get(), static_cast<double>(i),
                                         static_cast<double>(i) + 1.));
            });
            fwdpp::extensions::discrete_rec_model drm(
                recrate_region, std::move(f), std::move(w));
            vdrm_.emplace_back(std::move(drm));
        }
    return vdrm_;
}
