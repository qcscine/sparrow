/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef SPARROW_EXCEPTIONS_H
#define SPARROW_EXCEPTIONS_H

namespace Scine {
namespace Core {
class Calculator;
} // namespace Core
namespace Sparrow {

/** Time dependent exceptions */
class InvalidCalculatorType : public std::exception {
 public:
  const char* what() const noexcept final {
    return error.c_str();
  }

 protected:
  std::string error;
};

/**
 * @brief Exception thrown if a non-NDDO calculator is set.
 */
class InvalidCalculatorTypeForCIS : public InvalidCalculatorType {
 public:
  explicit InvalidCalculatorTypeForCIS(std::shared_ptr<Core::Calculator> method);
};

/**
 * @brief Exception thrown if a non-DFTB calculator is set.
 */
class InvalidCalculatorTypeForTDDFTB : public InvalidCalculatorType {
 public:
  explicit InvalidCalculatorTypeForTDDFTB(std::shared_ptr<Core::Calculator> method);
};

/**
 * @brief Exception for the case in which an unrestricted or triplet TD-DFTB calculation is asked for,
 *        but no spin constants are available.
 */
class SpinConstantsNotAvailableException : public std::exception {
  const char* what() const noexcept final {
    return "Triplet symmetry or Unrestricted TD-DFTB calculation requested, but some spin constants are unavailable.";
  }
};

/**
 * @brief Exception for invalid reference calculations.
 */
class InvalidReferenceCalculationException : public std::exception {
  const char* what() const noexcept final {
    return "Result class of reference calculation does not contain any energy. Remember running "
           "referenceCalculation().";
  }
};

/**
 * @brief Exception thrown if no reference calculator was set.
 */
class MissingReferenceCalculatorException : public std::exception {
  const char* what() const noexcept final {
    return "No reference calculator was set.";
  }
};

class InvalidSpinMultiplicityException : public std::exception {
  const char* what() const noexcept final {
    return "Invalid spin symmetry in excited states calculation from RHF reference.";
  }
};
} // namespace Sparrow
} // namespace Scine

#endif // SPARROW_EXCEPTIONS_H
