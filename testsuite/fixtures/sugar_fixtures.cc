#include "sugar_fixtures.hpp"
#include <fwdpp/poisson_xover.hpp>
#include <fwdpp/recbinder.hpp>

// MULTILOCUS
const int mlocuspop_popgenmut_fixture::nloci = 4;
const double mlocuspop_popgenmut_fixture::per_locus_rate = 5e-3;

std::vector<std::pair<double, double>>
mlocuspop_popgenmut_fixture::make_boundaries()
{
    std::vector<std::pair<double, double>> rv;
    for (int i = 0; i < nloci; ++i)
        {
            rv.emplace_back(i, i + 1);
        }
    return rv;
}

std::vector<mlocuspop_popgenmut_fixture::mutmodel>
mlocuspop_popgenmut_fixture::make_mutmodels()
{
    std::vector<mlocuspop_popgenmut_fixture::mutmodel> rv;
    for (int i = 0; i < nloci; ++i)
        {
            auto mm =
                [this,
                 i](fwdpp::flagged_mutation_queue &recbin,
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
                fwdpp::poisson_xover(per_locus_rate, i, i + 1),
                this->rng.get()));
        }
    return rv;
}

std::vector<std::function<std::vector<double>(void)>>
mlocuspop_popgenmut_fixture::fill_bound_recmodels()
{
    std::vector<std::function<std::vector<double>(void)>> rv;
    for (unsigned i = 0; i < nloci; ++i)
        {
            bound_recmodels.emplace_back(
                [this, i]() { return vdrm.at(i)(rng.get()); });
        }
    return rv;
}

std::vector<fwdpp::extensions::discrete_mut_model<
    mlocuspop_popgenmut_fixture::poptype::mcont_t>>
mlocuspop_popgenmut_fixture::fill_vdmm()
{
    // create a vector of extensions::discrete_mut_model
    std::vector<fwdpp::extensions::discrete_mut_model<poptype::mcont_t>> vdmm_;
    for (auto m : mutmodels)
        {
            vdmm_.emplace_back(
                fwdpp::extensions::discrete_mut_model<poptype::mcont_t>(
                    std::move(m), std::vector<double>(1, 1.0)));
        }

    return vdmm_;
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
