"""Facets of the X-ray polytope of doubly stochastic matrices: exact computations.

The polytope hat-P(k) is the set of anti-diagonal X-rays of k x k doubly stochastic matrices. The package contains the
certified facet test for the potential inequalities I(u), the facet criterion and the fiber of a facet, complete
enumerations of the rays of the normal fan (a Minkowski sum of simplices), and the checks of the statements about
step potentials, total support of staircases, padding, refinement, the non-step families and the height of facets.

Modules (see README.md and CLAIM_EVIDENCE_MAP.md for the claim each one supports):
  facet_test, facet_criterion, exact_graphs                  certified facet test, facet criterion, fibers
  fan_walk, vertex_walk, vertex_walk_cdd, small_k_fan        enumeration of the rays of the fan
  staircases, level_graph, hall_sets, boundary_classes       step potentials and total support of staircases
  cell_enumeration, sign_certificates                        statement E-flat (cells and certificate streams)
  attaining_structures, s1_certificates                      statement (S1) (catalogues and certificate trees)
  staircase_subclasses                                       classes of tails on which E-flat holds
  padding, refinement, block_lifts                           padding, cores, grid refinement, block lifts
  nonstep_families, potential_families                       non-step families, quadratic height, cubical fibers
"""
