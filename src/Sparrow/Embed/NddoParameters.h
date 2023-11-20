/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef INCLUDE_SPARROW_EMBED_NDDO_PARAMETERS_H
#define INCLUDE_SPARROW_EMBED_NDDO_PARAMETERS_H

#include "boost/filesystem.hpp"

//! Writes a cpp file containing Nddo parameters
void nddoParameters(boost::filesystem::path filename);

#endif
