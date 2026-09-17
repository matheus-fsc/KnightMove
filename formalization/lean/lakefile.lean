import Lake
open Lake DSL

package «knight_tour_lean» where
  leanOptions := #[
    ⟨`autoImplicit, false⟩
  ]

require mathlib from git
  "https://github.com/leanprover-community/mathlib4" @ "v4.16.0"

@[default_target]
lean_lib «KnightTour» where
  roots := #[`KnightTour]
