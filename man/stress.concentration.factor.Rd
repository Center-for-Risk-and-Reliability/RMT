\name{stress.concentration.factor}
\alias{stress.concentration.factor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Stress Concentration Factor Calculator
}
\description{
Computes the stress concentration factor \eqn{K_t} for a specified geometry and set of dimensions.
}
\usage{
stress.concentration.factor(dimensions, geometry, loadtype)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dimensions}{
Standard dimensions of the sample geometry based on several of the following: Width \code{W}, thickness \code{t}, semicircular notch radius \code{r}, U-notch depth \code{d}, and/or long diameter \code{d}.  These dimensions must be entered as a list.
}
  \item{geometry}{
Definition of the sample geometry from which a stress concentration factor is sought from.  The current list of geometries is as follows:

\strong{Rectangular Bar or Plate options}
  \itemize{
      \item \strong{RB1. Rectangular Bar with a Semi-circle Edge Notch} - Enter as \code{"rect_1semicirc_edge"} with dimensions width \code{W} and radius \code{r}.  Available \code{loadtype}: Tension (\code{1}) and Bending (\code{2}). (\emph{CONSTRAINT}: \eqn{r/W < 0.3})

      \figure{RB1.png}{options: alt="RB1 - Rectangular Bar with a Semi-circle Edge Notch"}

      \item \strong{RB2. Rectangular Bar with Opposite Semi-circle Edge Notches} - Enter as \code{"rect_2semicirc_edge"} with dimensions width \code{W} and radius \code{r}.  Available \code{loadtype}: Tension (\code{1}) and Bending (\code{2}). (\emph{CONSTRAINT}: \eqn{r/W < 0.5})

      \figure{RB2.png}{options: alt="RB2 - Rectangular Bar with Opposite Semi-circle Edge Notches"}

      \item \strong{RB3. Rectangular Bar with a U-Notch} - Enter as \code{"rect_1U_edge"} with dimensions width \code{W}, radius \code{r}, and U-depth \code{d}.  Available \code{loadtype}: Tension (\code{1}) and Bending (\code{2}).  (\emph{CONSTRAINT}: \eqn{0.5 < d/r < 20})

      \figure{RB3.png}{options: alt="RB3 - Rectangular Bar with a U-Notch"}

      \item \strong{RB4. Rectangular Bar with Opposite Edge U-Notches} - Enter as \code{"rect_2U_edge"} with dimensions width \code{W}, radius \code{r}, and U-depth \code{d}.   Available \code{loadtype}: Tension (\code{1}) and Bending (\code{2}).   (\emph{CONSTRAINT}: \eqn{0.02 \leq d/w \leq 0.25} and \eqn{0.1 \leq d/r \leq 50})

      \figure{RB4.png}{options: alt="RB4 - Rectangular Bar with Opposite Edge U-Notches"}

      \item \strong{RB6. Rectangular Bar with Opposite Edge V-Notches} - Enter as \code{"rect_2V_edge"} with dimensions width \code{W}, radius \code{r}, V-depth \code{d}, and V-angle \code{alp}.  Available \code{loadtype}: Tension (\code{1}).  (\emph{CONSTRAINT}: When \eqn{2d/W = 0.398}, \eqn{90^\circ \leq \alpha \leq 150^\circ} and When \eqn{2d/W = 0.667}, \eqn{60^\circ \leq \alpha \leq 150^\circ})

      \figure{RB6.png}{options: alt="RB6 - Rectangular Bar with Opposite Edge V-Notches"}

      \item \strong{RB7. Plate with a Circular Symmetric Hole} - Enter as \code{"rect_circ_symm_hole"} with dimensions width \code{W} and diameter \code{d}.    Available \code{loadtype}: Tension (\code{1}).  (\emph{CONSTRAINT}: \eqn{d \leq 0.9W})

      \figure{RB7.png}{options: alt="RB7 - Plate with a Circular Symmetric Hole"}

      \item \strong{RB8. Plate with a Elliptical Symmetric Hole} - Enter as \code{"rect_ellips_symm_hole"} with dimensions width \code{W}, short diameter \code{a}, and long diameter \code{d}.  Available \code{loadtype}: Tension (\code{1}).  (\emph{CONSTRAINT}: \eqn{1 < d/a < 8})

      \figure{RB8.png}{options: alt="RB8 - Plate with a Elliptical Symmetric Hole"}

      \item \strong{RB10. Plate with a Circular Offset Hole} - Enter as \code{"plate_circ_offset_hole"} with dimensions width \code{W}, radius \code{r}, and edge-to-hole center distance  \code{a}.  Available \code{loadtype}: Tension (\code{1}).  (\emph{CONSTRAINT}: \eqn{0 < r/a < 0.8})

      \figure{RB10.png}{options: alt="RB10 - Plate with a Circular Offset Hole"}

      \item \strong{RB13. Rectangular Bar with Infinite Rows of Semicircular Notches} - Enter as \code{"rect_inf_semicirc_edge"} with dimensions width \code{W}, radius \code{r}, and semi-circle spacing \code{b}.  Available \code{loadtype}: Tension (\code{1}) and Bending (\code{2}).  (\emph{CONSTRAINT}: \eqn{2r/W \leq 0.4} and \eqn{0 \leq 2r/b \leq 1})

      \figure{RB13.png}{options: alt="RB13 - Rectangular Bar with Infinite Rows of Semicircular Notches"}

      \item \strong{RB15. Rectangular Bar with Double Fillets} - Enter as \code{"rect_double_fillet"} with dimensions width \code{W}, fillet spacing \code{L}, inner width \code{d}, and radius \code{r}.  Available \code{loadtype}: Tension (\code{1}) and Bending (\code{2}).  (\emph{CONSTRAINT}: \eqn{0.1 \leq (W-d)/2r \leq 20.0})

      \figure{RB15.png}{options: alt="RB15 - Rectangular Bar with Double Fillets"}
  }

  \strong{Round or Shaft Options}
  \itemize{
      \item \strong{RS1. Round Shaft with a Single Fillet} - Enter as \code{"shaft_fillet_single"} with dimensions outer diameter \code{D}, inner diameter \code{d}, and radius \code{r}.  Available \code{loadtype}: Tension (\code{1}), Bending (\code{2}), and Torsion (\code{3}).  (\emph{CONSTRAINT}: \eqn{0.1 < \left(D - d\right)/2r < 20} for tension and bending, and \eqn{0.25 < \left(D - d\right)/2r < 4} for torsion)

      \figure{RS1.png}{options: alt="RS1 - Round Shaft with a Single Fillet"}

      \item \strong{RS4. Round Shaft with a U-Shaped or Semicircular Groove} - Enter as \code{"shaft_reg_U_groove"} with dimensions outer diameter \code{D}, depth of groove \code{t}, and radius of groove \code{r}.  Available \code{loadtype}: Bending (\code{2}) and Torsion (\code{3}).  (\emph{CONSTRAINT}: \eqn{0.25 < t/r < 50})

      \figure{RS4.png}{options: alt="RS4 - Round Shaft with a U-Shaped or Semicircular Groove"}

      \item \strong{RS5. Round Shaft with a V-Shaped Groove} - Enter as \code{"shaft_reg_V_groove"} with dimensions outer diameter \code{D}, depth of groove \code{t}, and radius of groove \code{r}.  Available \code{loadtype}: Bending (\code{2}) and Torsion (\code{3}).  (\emph{CONSTRAINT}: \eqn{'\alpha < 125^\circ}, \eqn{0.25 < t/r < 50}, and when \eqn{\alpha > 90^\circ}, \eqn{\frac{r}{D-2d} < 0.01})

      \figure{RS5.png}{options: alt="RS5 - Round Shaft with a V-Shaped Groove"}

      \item \strong{RS8. Cylindrical Tube with a Circular Hole} - Enter as \code{"shaft_cyl_tube_circ_hole"} with dimensions outer radius \code{R}, inner radius \code{r}, and diameter if hole \code{d}.  Available \code{loadtype}: Tension (\code{1}) and Bending (\code{2}).  (\emph{CONSTRAINT}: \eqn{'\frac{2\left(R-r\right)}{R+r} < 125^\circ})

      \figure{RS8.png}{options: alt="RS8 - Cylindrical Tube with a Circular Hole"}

      \item \strong{RS9. Round Bar with as Shallow U-Shaped (dog-bone) Groove} - Enter as \code{"shaft_wideU_groove"} with dimensions outer diameter \code{D}, inner diameter \code{d}, and radius of groove \code{r}.  Available \code{loadtype}: Tension (\code{1}), Bending (\code{2}), and Torsion (\code{3}).  (\emph{CONSTRAINT}: \eqn{1.005 < D/d < 1.1} and \eqn{0.3 < r/d < 1})

      \figure{RS9.png}{options: alt="RS9 - Round Bar with as Shallow U-Shaped (dog-bone) Groove"}
  }

  \strong{Infinite/Semi-Infinite Plate Options}
  \itemize{
      \item \strong{ISI2. SInfinite Plate with Two Circular Holes} - Enter as \code{"inf_plate_two_circ_holes"} with dimensions distance between two holes \code{a} and radius \code{r}.  Available \code{loadtype}: Parallel Tension to Holes (\code{1}) and Perpendicular Tension to Holes (\code{2}).  (\emph{CONSTRAINT}: \eqn{0 \leq \frac{2r}{a} < 10})

      \figure{ISI2.png}{options: alt="ISI2 - Infinite Plate with Two Circular Holes"}

      \item \strong{ISI4. Semi-Infinite Plate with a Circular Offset Hole} - Enter as \code{"semi_inf_plate_circ_offset_hole"} with dimensions edge to hole center distance \code{a} and radius \code{r}.  Available \code{loadtype}: Tension (\code{1}).  (\emph{CONSTRAINT}: \eqn{r < a})

      \figure{ISI4.png}{options: alt="ISI4 - Semi-Infinite Plate with a Circular Offset Hole"}
  }

  \strong{Pressurized Round or Shaft Options}
  \itemize{
      \item \strong{PRS1. Pressurized (Internal) Cylindrical Tube with a Circular Hole} - Enter as \code{"shaft_cyl_tube_circ_hole_pressurized"} with dimensions outer radius \code{R}, inner radius \code{r}, and diameter if hole \code{d}.  Available \code{loadtype}: Tension (\code{1}) and Bending (\code{2}).  (\emph{CONSTRAINT}: \eqn{0 \leq \frac{\left(R+r\right)}{\left(R-r\right)} \leq 12.52068} for Membrane effects and \eqn{0 \leq \frac{\left(R+r\right)}{\left(R-r\right)} \leq 9.828731} for Membrane plus bending effects)

      \figure{PRS1.png}{options: alt="PRS1 - Pressurized (Internal) Cylindrical Tube with a Circular Hole"}

      \item \strong{PRS2. Pressurized Spherical Shell with Elliptical Hole} - Enter as \code{"shell_spherical_ellips_hole_pressurized"} with dimensions outer radius of sphere \code{R}, inner radius of sphere \code{r}, semi-major radius of ellipse \code{a}, and semi-minor radius of ellipse \code{b}.  Available \code{loadtype}: Tension (\code{1}).  (\emph{CONSTRAINT}: \eqn{r<R}, \eqn{a<r}, and \eqn{b<r})

      \figure{PRS2.png}{options: alt="PRS2 - PPressurized Spherical Shell with Elliptical Hole"}
  }

}
  \item{loadtype}{
The description of the load on the sample: Enter \code{1} for tension (default), enter \code{2} for bending, and enter \code{3} for traverse bending.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Computes the stress concentration factor \eqn{K_t}.
}
\references{
Pikley, W., Pikley, D., & Bi, Z. (2020). Peterson's Stress Concentration Factors (4th Edition). Wiley
}
\author{
Reuel Smith (2022).  Additional geometry code contributions by Camille Levine and Colin Schell.
}
\note{
Additional geometry types may be added to the current list in future versions.  Please see available \code{loadtype} options for current list.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
stress.concentration.factor(list(W=100,r=10),"rect_2semicirc_edge")

stress.concentration.factor(list(W=100,r=10),"rect_1semicirc_edge")

stress.concentration.factor(list(r = 20, W = 100),"rect_2semicirc_edge")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
