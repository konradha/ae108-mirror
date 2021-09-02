#include <complex>
#include <petscsys.h>
#include <type_traits>

int main() {
#ifdef AE108_PETSC_COMPLEX
  static_assert(std::is_same<PetscScalar, std::complex<double>>::value, "");
#endif

#ifdef AE108_PETSC_REAL
  static_assert(std::is_same<PetscScalar, double>::value, "");
#endif
}