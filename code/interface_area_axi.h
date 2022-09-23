// Modification of the standard interface_area from fractions.h to work for
// an axisymmetric simulation. As the standard version would calculate the 
// arc length for an axi simulation by using the Theorem of Pappus by 
// multiplying the arc length of the segment by distance travelled by the 
// centroid when rotating (which is 2*pi*r) and in this case r is the y 
// coordinate.

double interface_area_axi (scalar c)
 {
  double area = 0.;
  foreach (reduction(+:area))
    if (c[] > 1e-6 && c[] < 1. - 1e-6) {
      coord n = interface_normal (point, c), p;
      double alpha = plane_alpha (c[], n);
      area += pow(Delta, dimension - 1)*plane_area_center (n, alpha, &p) *2*M_PI*y;
    }
  return area;
}