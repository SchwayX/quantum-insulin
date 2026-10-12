"""Collect every number quoted in AUDIT.md / README.md from the raw result files (no new calculations,
apart from the clearly labelled derived estimates at the bottom).

    python audit/2026-10-12/scripts/key_numbers.py audit/2026-10-12/results
"""
import json
import os
import sys

EPS_CHEM = 1.6e-3  # Eh, ~1 kcal/mol ("chemical accuracy") used for the shot-count estimate


def load(d, name):
    with open(os.path.join(d, name)) as fh:
        return json.load(fh)


def main(d):
    h = load(d, "h2_audit.json")
    c = load(d, "cys_audit.json")
    f = load(d, "cys_followup.json")
    s = load(d, "scaling.json")
    p0, p1 = h["single_points"]
    K = {"source_files": ["h2_audit.json", "cys_audit.json", "cys_followup.json", "scaling.json"]}
    K["software"] = {"h2": h["software"], "cysteine": c["software"]}
    K["H2_R0.7414"] = {k: p0[k] for k in ("E_nuc", "E_RHF", "E_FCI_pyscf", "E_JW_exact_N2_Sz0", "E_CI2x2_analytic",
                                          "n_qubits", "n_pauli_terms_incl_identity", "pauli_1norm_excl_identity",
                                          "lowest_eigenvalue_by_particle_number")}
    K["H2_R0.735"] = {k: p1[k] for k in ("E_nuc", "E_RHF", "E_FCI_pyscf")}
    K["H2_checks"] = {"points": h["n_points_checked"], "worst_case": h["worst_case_checks"],
                      "all_routes_agree": h["all_routes_agree"], "variational_bound_ok": h["variational_bound_ok"],
                      "global_min_always_N2": h["global_min_always_in_N2"],
                      "R_min_RHF_A": h["minimum_RHF_A"], "R_min_FCI_A": h["minimum_FCI_A"],
                      "E_two_H_atoms": h["E_two_H_atoms_STO3G"],
                      "E_FCI_3.0A": h["scan"][-1]["E_FCI_pyscf"], "E_RHF_3.0A": h["scan"][-1]["E_RHF"],
                      "wall_time_s": h["wall_time_s"]}
    K["cysteine_counts"] = c["counts"]
    K["cysteine_conformers"] = [{"seed": x["seed(chi1,chi2,psi)"], "E_start": x["E_idealized_start"],
                                 "E_opt": x["E_optimized"], "opt_converged": x["opt_converged"],
                                 "max_grad": x["max_abs_gradient_Eh_per_bohr"], "opt_time_s": x["opt_time_s"],
                                 "bonds": x["geometry"]["n_bonds_detected"], "stability": x["stability"]}
                                for x in c["conformers"]]
    K["cysteine_energy"] = {"E_range_optimized": c["E_range_optimized"],
                            "spread_optimized_mEh": 1e3 * (c["E_range_optimized"][1] - c["E_range_optimized"][0]),
                            "E_range_idealized_starts": c["E_range_idealized_starts"],
                            "E_lowest": c["lowest_conformer_E"], "historical_minus_lowest": c["historical_minus_lowest"]}
    orb = [x["orbitals"] for x in c["conformers"]]
    K["cysteine_orbitals_range"] = {k: [min(o[k] for o in orb), max(o[k] for o in orb)]
                                    for k in ("E_HOMO_Eh", "E_LUMO_Eh", "gap_eV", "koopmans_IP_eV",
                                              "HOMO_mulliken_weight_on_S")}
    K["cysteine_frequencies_lowest"] = c["lowest_conformer_frequencies"]
    hyp = {}
    for k, v in c["error_hypotheses_vs_lowest_conformer"].items():
        if "dE" in v:
            hyp[k] = v["dE"]
        else:
            hyp[k] = {kk: {"min": vv["dE_min"], "mean": vv["dE_mean"], "max": vv["dE_max"]} for kk, vv in v.items()}
    K["error_hypotheses_dE_vs_lowest"] = hyp
    K["MP2_correlation"] = c["MP2_correlation_full_space"]
    K["active_spaces"] = c["active_spaces"]
    K["full_space_one_norm"] = c["full_space_one_norm"]
    K["external_instability"] = f["external_stability"]
    K["pauli_terms"] = {"by_norb": f["pauli_terms_by_norb"], "fit": f["poly4_fit_on_norb_2_to_9"],
                        "extrapolated": f["extrapolated_pauli_terms"]}
    K["scaling"] = s["systems"]
    # ---- derived estimates (assumption-laden; stated as such in the documents)
    lam = {"H2_R0.7414": p0["pauli_1norm_excl_identity"], "cysteine_CAS(12,12)": c["active_spaces"][-1]["one_norm_koridon"],
           "cysteine_full_46_orbitals": c["full_space_one_norm"]}
    K["derived"] = {
        "shots_one_norm_bound_(lambda/eps)^2_eps=1.6mEh": {k: (v / EPS_CHEM) ** 2 for k, v in lam.items()},
        "statevector_complex128_GiB": {n: 2 ** n * 16 / 2 ** 30 for n in (20, 24, 26, 28, 30, 32)},
        "density_matrix_complex128_GiB": {n: 4 ** n * 16 / 2 ** 30 for n in (10, 12, 13, 14)},
        "CAS12_corr_fraction_of_frozen_core_MP2": c["active_spaces"][-1]["E_corr_captured"] / c["MP2_correlation_full_space"]["frozen_core"],
    }
    out = os.path.join(d, "key_numbers.json")
    with open(out, "w") as fh:
        json.dump(K, fh, indent=1)
    print("wrote", out)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
