// _LICENSE_HEADER_
//
// Copyright (C) 2019 - 2025.
// Terms register on the GPL-3.0 license.
//
// This file can be redistributed and/or modified under the license terms.
//
// See top level LICENSE file for more details.
//
// This file can be used citing references in CITATION.cff file.

#include "MeshUtilities.hpp"

namespace Gedim
{
  // ***************************************************************************
  inline constexpr std::array<std::array<int,16>,256> TriTable =
  {
    {
      {{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
      {{0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
      {{0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
      {{1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
      {{1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
      {{0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
      {{9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
      {{2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1}},
      {{3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
      {{0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
      {{1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
      {{1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1}},
      {{3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
      {{0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1}},
      {{3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1}},
      {{9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
      {{4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
      {{4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
      {{0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
      {{4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1}},
      {{1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
      {{3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1}},
      {{9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1}},
      {{2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1}},
      {{8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
      {{11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1}},
      {{9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1}},
      {{4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1}},
      {{3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1}},
      {{1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1}},
      {{4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1}},
      {{4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1}},
      {{9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
      {{9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
      {{0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
      {{8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1}},
      {{1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
      {{3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1}},
      {{5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1}},
      {{2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1}},
      {{9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
      {{0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1}},
      {{0, 5, 4, 0, 1, 5, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1}},
      {{2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1}},
      {{10, 3, 11, 10, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1}},
      {{4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1, -1, -1, -1}},
      {{5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1}},
      {{5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1}},
      {{9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
      {{9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1}},
      {{0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1}},
      {{1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
      {{9, 7, 8, 9, 5, 7, 10, 1, 2, -1, -1, -1, -1, -1, -1, -1}},
      {{10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1, -1}},
      {{8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1}},
      {{2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1}},
      {{7, 9, 5, 7, 8, 9, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1}},
      {{9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1}},
      {{2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1}},
      {{11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1}},
      {{9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1, -1, -1, -1}},
      {{5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1}},
      {{11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1}},
      {{11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
      {{10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
      {{0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
      {{9, 0, 1, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
      {{1, 8, 3, 1, 9, 8, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1}},
      {{1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
      {{1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1, -1}},
      {{9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1}},
      {{5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1}},
      {{2, 3, 11, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
      {{11, 0, 8, 11, 2, 0, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1}},
      {{0, 1, 9, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1}},
      {{5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1, -1, -1, -1}},
      {{6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1}},
      {{0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1}},
      {{3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1}},
      {{6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1}},
      {{5, 10, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
      {{4, 3, 0, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1}},
      {{1, 9, 0, 5, 10, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1}},
      {{10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1, -1}},
      {{6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1}},
      {{1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1, -1}},
      {{8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1, -1}},
      {{7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1}},
      {{3, 11, 2, 7, 8, 4, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1}},
      {{5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1, -1, -1, -1}},
      {{0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1}},
      {{9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1}},
      {{8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1, -1, -1, -1}},
      {{5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1}},
      {{0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1}},
      {{6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1, -1, -1, -1}},
      {{10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
      {{4, 10, 6, 4, 9, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1}},
      {{10, 0, 1, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1}},
      {{8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1}},
      {{1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1}},
      {{3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1, -1, -1, -1}},
      {{0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
      {{8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1}},
      {{10, 4, 9, 10, 6, 4, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1}},
      {{0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1, -1, -1, -1}},
      {{3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1, -1, -1, -1}},
      {{6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1}},
      {{9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1, -1, -1, -1}},
      {{8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1}},
      {{3, 11, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1, -1}},
      {{6, 4, 8, 11, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
      {{7, 10, 6, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1}},
      {{0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1, -1, -1, -1}},
      {{10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1, -1}},
      {{10, 6, 7, 10, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1}},
      {{1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1, -1}},
      {{2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1}},
      {{7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1}},
      {{7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
      {{2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1, -1, -1, -1}},
      {{2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1}},
      {{1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1}},
      {{11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1, -1, -1, -1}},
      {{8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1}},
      {{0, 9, 1, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
      {{7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1, -1, -1, -1}},
      {{7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
      {{7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
      {{3, 0, 8, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
      {{0, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
      {{8, 1, 9, 8, 3, 1, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1}},
      {{10, 1, 2, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
      {{1, 2, 10, 3, 0, 8, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1}},
      {{2, 9, 0, 2, 10, 9, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1}},
      {{6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, -1, -1, -1, -1}},
      {{7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
      {{7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1, -1, -1, -1, -1}},
      {{2, 7, 6, 2, 3, 7, 0, 1, 9, -1, -1, -1, -1, -1, -1, -1}},
      {{1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1, -1, -1, -1}},
      {{10, 7, 6, 10, 1, 7, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1}},
      {{10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1, -1, -1, -1}},
      {{0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1, -1, -1, -1}},
      {{7, 6, 10, 7, 10, 8, 8, 10, 9, -1, -1, -1, -1, -1, -1, -1}},
      {{6, 8, 4, 11, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
      {{3, 6, 11, 3, 0, 6, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1}},
      {{8, 6, 11, 8, 4, 6, 9, 0, 1, -1, -1, -1, -1, -1, -1, -1}},
      {{9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1, -1, -1, -1}},
      {{6, 8, 4, 6, 11, 8, 2, 10, 1, -1, -1, -1, -1, -1, -1, -1}},
      {{1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, -1, -1, -1, -1}},
      {{4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, -1, -1, -1, -1}},
      {{10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, -1}},
      {{8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1}},
      {{0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
      {{1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1, -1, -1, -1}},
      {{1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1}},
      {{8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1, -1, -1, -1}},
      {{10, 1, 0, 10, 0, 6, 6, 0, 4, -1, -1, -1, -1, -1, -1, -1}},
      {{4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, -1}},
      {{10, 9, 4, 6, 10, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
      {{4, 9, 5, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
      {{0, 8, 3, 4, 9, 5, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1}},
      {{5, 0, 1, 5, 4, 0, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1}},
      {{11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1, -1, -1, -1}},
      {{9, 5, 4, 10, 1, 2, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1}},
      {{6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, -1, -1, -1, -1}},
      {{7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, -1, -1, -1, -1}},
      {{3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, -1}},
      {{7, 2, 3, 7, 6, 2, 5, 4, 9, -1, -1, -1, -1, -1, -1, -1}},
      {{9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1, -1, -1, -1}},
      {{3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1, -1, -1, -1}},
      {{6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, -1}},
      {{9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, -1, -1, -1, -1}},
      {{1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, -1}},
      {{4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, -1}},
      {{7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, -1, -1, -1, -1}},
      {{6, 9, 5, 6, 11, 9, 11, 8, 9, -1, -1, -1, -1, -1, -1, -1}},
      {{3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1, -1, -1, -1}},
      {{0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1, -1, -1, -1}},
      {{6, 11, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1, -1, -1, -1, -1}},
      {{1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, -1, -1, -1, -1}},
      {{0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, -1}},
      {{11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, -1}},
      {{6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, -1, -1, -1, -1}},
      {{5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1, -1, -1, -1}},
      {{9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1, -1, -1, -1, -1}},
      {{1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, -1}},
      {{1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
      {{1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, -1}},
      {{10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, -1, -1, -1, -1}},
      {{0, 3, 8, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
      {{10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
      {{11, 5, 10, 7, 5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
      {{11, 5, 10, 11, 7, 5, 8, 3, 0, -1, -1, -1, -1, -1, -1, -1}},
      {{5, 11, 7, 5, 10, 11, 1, 9, 0, -1, -1, -1, -1, -1, -1, -1}},
      {{10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, -1, -1, -1, -1}},
      {{11, 1, 2, 11, 7, 1, 7, 5, 1, -1, -1, -1, -1, -1, -1, -1}},
      {{0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, -1, -1, -1, -1}},
      {{9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1, -1, -1, -1}},
      {{7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, -1}},
      {{2, 5, 10, 2, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1}},
      {{8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1, -1, -1, -1}},
      {{9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, -1, -1, -1, -1}},
      {{9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, -1}},
      {{1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
      {{0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1, -1, -1, -1, -1}},
      {{9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1, -1, -1, -1, -1}},
      {{9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
      {{5, 8, 4, 5, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1}},
      {{5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1, -1, -1, -1}},
      {{0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, -1, -1, -1, -1}},
      {{10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, -1}},
      {{2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1, -1, -1, -1}},
      {{0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, -1}},
      {{0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, -1}},
      {{9, 4, 5, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
      {{2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1, -1, -1, -1}},
      {{5, 10, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1, -1, -1, -1, -1}},
      {{3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, -1}},
      {{5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1, -1, -1, -1}},
      {{8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1, -1, -1, -1, -1}},
      {{0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
      {{8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1, -1, -1, -1}},
      {{9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
      {{4, 11, 7, 4, 9, 11, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1}},
      {{0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, -1, -1, -1, -1}},
      {{1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1, -1, -1, -1}},
      {{3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, -1}},
      {{4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1, -1, -1, -1}},
      {{9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, -1}},
      {{11, 7, 4, 11, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1, -1}},
      {{11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, -1, -1, -1, -1}},
      {{2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1, -1, -1, -1}},
      {{9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, -1}},
      {{3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, -1}},
      {{1, 10, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
      {{4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1, -1, -1, -1, -1}},
      {{4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1, -1, -1, -1}},
      {{4, 0, 3, 7, 4, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
      {{4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
      {{9, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
      {{3, 0, 9, 3, 9, 11, 11, 9, 10, -1, -1, -1, -1, -1, -1, -1}},
      {{0, 1, 10, 0, 10, 8, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1}},
      {{3, 1, 10, 11, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
      {{1, 2, 11, 1, 11, 9, 9, 11, 8, -1, -1, -1, -1, -1, -1, -1}},
      {{3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, -1, -1, -1, -1}},
      {{0, 2, 11, 8, 0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
      {{3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
      {{2, 3, 8, 2, 8, 10, 10, 8, 9, -1, -1, -1, -1, -1, -1, -1}},
      {{9, 10, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
      {{2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, -1, -1, -1, -1}},
      {{1, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
      {{1, 3, 8, 9, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
      {{0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
      {{0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
      {{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
    }
  };
  // ***************************************************************************
  using Triangle = std::array<unsigned int, 3>;
  inline constexpr unsigned int None = std::numeric_limits<unsigned int>::max();
  inline constexpr std::array<std::array<int, 2>,12> EdgeEnds = {{{{0,1}},{{1,2}},{{2,3}},{{3,0}},{{4,5}},{{5,6}},{{6,7}},{{7,4}},{{0,4}},{{1,5}},{{2,6}},{{3,7}}}};
  inline constexpr std::array<std::array<int, 3>,8> Bits =
  {
    {
      {{0,0,0}},{{1,0,0}},{{1,1,0}},{{0,1,0}},{{0,0,1}},{{1,0,1}},{{1,1,1}},{{0,1,1}}
    }
  };

  struct Cube
  {
      std::array<unsigned int, 8> vertices;
      std::array<unsigned int, 12> edges;
      Eigen::Vector3d lo, hi;
  };

  inline Cube BuildCube(const Gedim::GeometryUtilities& geometry_utilities,
                        const Gedim::IMeshDAO& mesh,
                        unsigned int cell)
  {
    if (mesh.Cell3DNumberVertices(cell) != 8 ||
        mesh.Cell3DNumberEdges(cell)!=12)
      throw std::runtime_error("Expected hexahedron with 8 vertices and 12 edges");

    Cube c; c.vertices.fill(None); c.edges.fill(None);
    c.lo.setConstant(std::numeric_limits<double>::infinity());
    c.hi.setConstant(-std::numeric_limits<double>::infinity());
    for (unsigned int k = 0; k < 8; ++k)
    {
      const auto coordinate = mesh.Cell0DCoordinates(mesh.Cell3DVertex(cell, k));
      for(unsigned int a=0; a < 3; ++a)
      {
        c.lo[a] = std::min(c.lo[a], coordinate[a]);
        c.hi[a] = std::max(c.hi[a], coordinate[a]);
      }
    }

    const Eigen::Vector3d extent = c.hi - c.lo;
    if (geometry_utilities.IsValueLowerOrEqual(extent.x(),
                                               0.0,
                                               geometry_utilities.Tolerance1D()) ||
        geometry_utilities.IsValueLowerOrEqual(extent.y(),
                                               0.0,
                                               geometry_utilities.Tolerance1D()) ||
        geometry_utilities.IsValueLowerOrEqual(extent.z(),
                                               0.0,
                                               geometry_utilities.Tolerance1D()))
      throw std::runtime_error("Zero cube extent");

    for (unsigned int k = 0; k < 8; ++k)
    {
      const unsigned int cell0D_id = mesh.Cell3DVertex(cell, k);
      const auto coordinate = mesh.Cell0DCoordinates(cell0D_id);
      int bits[3];

      for (unsigned int a = 0; a < 3; ++a)
      {
        const double t = (coordinate[a] - c.lo[a]) / extent[a];

        if (geometry_utilities.IsValueZero(t,
                                           geometry_utilities.Tolerance1D()))
          bits[a] = 0;
        else if (geometry_utilities.IsValueZero((t - 1.0),
                                                geometry_utilities.Tolerance1D()))
          bits[a] = 1;
        else
          throw std::runtime_error("Cell is not axis-aligned Cartesian");
      }

      int local = -1;
      for (unsigned int j = 0; j < 8; ++j)
      {
        if (Bits[j][0] == bits[0] &&
            Bits[j][1] == bits[1] &&
            Bits[j][2] == bits[2])
        {
          local=j;
          break;
        }
      }

      if (local < 0 ||
          c.vertices[local] != None)
        throw std::runtime_error("Invalid cube vertex ordering");

      c.vertices[local] = cell0D_id;
    }

    for (unsigned int e = 0; e< 12; ++e)
    {
      const unsigned int a = c.vertices[EdgeEnds[e][0]];
      const unsigned int b = c.vertices[EdgeEnds[e][1]];

      for (unsigned int j = 0; j < 12; ++j)
      {
        const unsigned int cell1D_id = mesh.Cell3DEdge(cell,j);
        const unsigned int origin_id = mesh.Cell1DOrigin(cell1D_id);
        const unsigned int end_id = mesh.Cell1DEnd(cell1D_id);

        if ((origin_id==a && end_id==b) ||
            (origin_id==b && end_id==a))
        {
          c.edges[e] = cell1D_id;
          break;
        }
      }

      if (c.edges[e]==None)
        throw std::runtime_error("Missing cube edge");
    }

    return c;
  }
  // ***************************************************************************
  struct TriangleHash
  {
      std::size_t operator()(const Triangle& a) const noexcept
      {
        std::size_t h=0;
        for(auto x:a) h^=std::hash<unsigned int>{}(x)+0x9e3779b9+(h<<6)+(h>>2);
        return h;
      }
  };
  // ***************************************************************************
  inline Eigen::Vector3d GradientAt(const Cube& cube,
                                    const std::array<double, 8>& faces,
                                    const Eigen::Vector3d& point)
  {
    const Eigen::Vector3d t = ((point - cube.lo).array() /
                               (cube.hi - cube.lo).array()).matrix();
    Eigen::Vector3d g = Eigen::Vector3d::Zero();
    for (unsigned int i = 0; i < 8; ++i)
    {
      for (unsigned int a = 0; a < 3; ++a)
      {
        double w = (Bits[i][a] ? 1.0 : -1.0) / (cube.hi[a] - cube.lo[a]);
        for (unsigned int b = 0; b < 3; ++b)
        {
          if (b != a)
            w *= Bits[i][b] ? t[b] : (1.0 - t[b]);
        }

        g[a] += faces[i] * w;
      }
    }
    return g;
  }
  // ***************************************************************************
  inline std::vector<Triangle> CreateTriangles(const Gedim::GeometryUtilities& geometry_utilities,
                                               const Eigen::VectorXd& cell0Ds_surface_value,
                                               const std::vector<Gedim::MeshUtilities::MarchingCubesResult::SurfacePositions>& cell0Ds_surface_position,
                                               const std::vector<Gedim::MeshUtilities::IntersectMeshWithSurfaceResult::SurfacePointData>& surface_intersection_points,
                                               const std::vector<unsigned int>& cell0Ds_surface_points,
                                               const std::vector<unsigned int>& cell1Ds_surface_points,
                                               const Gedim::IMeshDAO& mesh)
  {
    using Position = Gedim::MeshUtilities::MarchingCubesResult::SurfacePositions;

    if (static_cast<unsigned int>(cell0Ds_surface_value.size()) != mesh.Cell0DTotalNumber() ||
        cell0Ds_surface_position.size() != mesh.Cell0DTotalNumber() ||
        cell0Ds_surface_points.size() != mesh.Cell0DTotalNumber() ||
        cell1Ds_surface_points.size() != mesh.Cell1DTotalNumber())
      throw std::invalid_argument("Surface arrays do not match mesh dimensions");

    std::vector<Triangle> output;
    std::unordered_set<Triangle,TriangleHash> seen;
    for (unsigned int cell = 0; cell < mesh.Cell3DTotalNumber(); ++cell)
    {
      const Cube c = BuildCube(geometry_utilities,
                               mesh,
                               cell);
      unsigned mask=0;
      std::array<double,8> values;
      for (unsigned int i = 0; i < 8; ++i)
      {
        const unsigned int v = c.vertices[i];
        values[i] = cell0Ds_surface_value[v];
        if (cell0Ds_surface_position[v] == Position::Negative)
          mask |= 1u<<i;
      }

      // Symbolic convention: OnSurface belongs to the nonnegative side.
      if(mask==0 || mask==255)
        continue;

      const auto& row = TriTable[mask];
      for (unsigned int j = 0; j < 16 && row[j] >=0; j+=3)
      {
        if (j + 2 >= 16 ||
            row[j + 1] < 0 || row[j+2] < 0)
          throw std::runtime_error("Malformed triangle table");

        Triangle tri{};
        for(unsigned int k = 0; k < 3; ++k)
        {
          int e = row[j + k];
          unsigned int a = c.vertices[EdgeEnds[e][0]];
          unsigned int b = c.vertices[EdgeEnds[e][1]];
          const bool az = (cell0Ds_surface_position[a] == Position::OnSurface);
          const bool bz = (cell0Ds_surface_position[b] == Position::OnSurface);

          // If both endpoints are OnSurface, use a deterministic endpoint.
          // Such a case is not fully defined by the classic table.
          unsigned int id = (az && bz) ?
                              cell0Ds_surface_points[std::min(a,b)] :
                            az ? cell0Ds_surface_points[a] :
                                 bz ? cell0Ds_surface_points[b] : cell1Ds_surface_points[c.edges[e]];
          if (id == None ||
              id >= surface_intersection_points.size())
            throw std::runtime_error("Required precomputed surface point is missing");

          tri[k] = id;
        }

        if (tri[0] == tri[1] ||
            tri[1] == tri[2] ||
            tri[0] == tri[2])
          continue;

        const auto& p0 = surface_intersection_points[tri[0]].coordinate;
        const auto& p1 = surface_intersection_points[tri[1]].coordinate;
        const auto& p2 = surface_intersection_points[tri[2]].coordinate;
        if (p0.size() < 3 ||
            p1.size() < 3 ||
            p2.size() < 3)
          throw std::runtime_error("Expected 3D intersection points");

        const Eigen::Vector3d a = p0.head<3>(), b = p1.head<3>(), d = p2.head<3>();
        const Eigen::Vector3d n = (b - a).cross(d - a);
        const double scale2 = std::max(
                                {
                                  (b - a).squaredNorm(),
                                  (d - a).squaredNorm(),
                                  (a - d).squaredNorm()
                                });
        if(geometry_utilities.IsValueZero(scale2,
                                          geometry_utilities.Tolerance1DSquared()) ||
           n.norm() <= geometry_utilities.Tolerance2D() * scale2)
          continue;

        const Eigen::Vector3d gradient = GradientAt(c,
                                                    values,
                                                    (a + b + d) / 3.0);

        if (geometry_utilities.IsValuePositive(gradient.squaredNorm(),
                                               geometry_utilities.Tolerance1DSquared())
            &&
            geometry_utilities.IsValueNegative(n.dot(gradient),
                                               geometry_utilities.Tolerance1D()))
          std::swap(tri[1], tri[2]);
        Triangle key = tri;
        std::sort(key.begin(),
                  key.end());

        if(seen.insert(key).second)
          output.push_back(tri);
      }
    }

    return output;
  }
  // ***************************************************************************
  Gedim::MeshUtilities::IntersectMeshWithSurfaceResult Gedim::MeshUtilities::IntersectMeshWithSurface(const Gedim::GeometryUtilities &geometry_utilities,
                                                                                                      const Eigen::VectorXd &cell0Ds_surface_value,
                                                                                                      const std::vector<Gedim::MeshUtilities::MarchingCubesResult::SurfacePositions> &cell0Ds_surface_position,
                                                                                                      const Gedim::IMeshDAO &mesh)
  {
    using SurfacePosition = Gedim::MeshUtilities::MarchingCubesResult::SurfacePositions;

    Gedim::MeshUtilities::IntersectMeshWithSurfaceResult result;

    result.cell0Ds_surface_points.resize(mesh.Cell0DTotalNumber(),
                                         Gedim::MeshUtilities::IntersectMeshWithSurfaceResult::None);
    result.cell1Ds_surface_points.resize(mesh.Cell1DTotalNumber(),
                                         Gedim::MeshUtilities::IntersectMeshWithSurfaceResult::None);
    result.surface_intersection_points.reserve(mesh.Cell1DTotalNumber() +
                                               mesh.Cell0DTotalNumber());

    for (unsigned int e = 0; e < mesh.Cell1DTotalNumber(); ++e)
    {
      const auto origin_id = mesh.Cell1DOrigin(e);
      const auto end_id = mesh.Cell1DEnd(e);

      const auto origin_coordinate = mesh.Cell1DOriginCoordinates(e);
      const auto end_coordinate = mesh.Cell1DEndCoordinates(e);

      const auto origin_surface_position = cell0Ds_surface_position.at(origin_id);
      const auto end_surface_position = cell0Ds_surface_position.at(end_id);

      switch (origin_surface_position)
      {
        case SurfacePosition::OnSurface:
        {
          if (result.cell0Ds_surface_points.at(origin_id) ==
              Gedim::MeshUtilities::IntersectMeshWithSurfaceResult::None)
          {
            result.cell0Ds_surface_points.at(origin_id) =
                result.surface_intersection_points.size();
            result.surface_intersection_points.push_back({
                                                           origin_coordinate,
                                                           Gedim::MeshUtilities::IntersectMeshWithSurfaceResult::SurfacePointData::Types::Vertex,
                                                           origin_id
                                                         });
          }
        }
          break;
        case SurfacePosition::Positive:
        case SurfacePosition::Negative:
          break;
        default:
          throw std::runtime_error("Position not managed");
      }

      switch (end_surface_position)
      {
        case SurfacePosition::OnSurface:
        {
          if (result.cell0Ds_surface_points.at(end_id) ==
              Gedim::MeshUtilities::IntersectMeshWithSurfaceResult::None)
          {
            result.cell0Ds_surface_points.at(end_id) =
                result.surface_intersection_points.size();
            result.surface_intersection_points.push_back({
                                                           end_coordinate,
                                                           Gedim::MeshUtilities::IntersectMeshWithSurfaceResult::SurfacePointData::Types::Vertex,
                                                           end_id
                                                         });
          }
        }
          break;
        case SurfacePosition::Positive:
        case SurfacePosition::Negative:
          break;
        default:
          throw std::runtime_error("Position not managed");
      }

      if (origin_surface_position != end_surface_position &&
          origin_surface_position != SurfacePosition::OnSurface &&
          end_surface_position != SurfacePosition::OnSurface)
      {
        const double origin_surface_value = cell0Ds_surface_value[origin_id];
        const double end_surface_value = cell0Ds_surface_value[end_id];

        const Eigen::Vector3d intersection = origin_coordinate +
                                             origin_surface_value *
                                             (end_coordinate - origin_coordinate) /
                                             (origin_surface_value - end_surface_value);

        result.cell1Ds_surface_points.at(e) =
            result.surface_intersection_points.size();
        result.surface_intersection_points.push_back({
                                                       intersection,
                                                       Gedim::MeshUtilities::IntersectMeshWithSurfaceResult::SurfacePointData::Types::Edge,
                                                       e
                                                     });

      }
    }

    result.surface_intersection_points.shrink_to_fit();

    return result;
  }
  // ***************************************************************************
  Gedim::MeshUtilities::MarchingCubesResult Gedim::MeshUtilities::MarchingCubes(
      const Gedim::GeometryUtilities &geometry_utilities,
      const std::function<Eigen::VectorXd(const Eigen::MatrixXd &)> &surface_implicit_equation,
      const Gedim::IMeshDAO &mesh)
  {
    using SurfacePosition = Gedim::MeshUtilities::MarchingCubesResult::SurfacePositions;

    Gedim::MeshUtilities::MarchingCubesResult result;

    result.cell0Ds_surface_value = surface_implicit_equation(mesh.Cell0DsCoordinates());

    {
      auto marking_function = [&geometry_utilities](const Eigen::VectorXd &points_surface_value) {
        std::vector<SurfacePosition> marks(points_surface_value.size());

        for (unsigned int p = 0; p < points_surface_value.size(); p++)
        {
          const auto surface_value = points_surface_value[p];

          if (geometry_utilities.IsValuePositive(surface_value, geometry_utilities.Tolerance1D()))
          {
            marks[p] = SurfacePosition::Positive;
          }
          else if (geometry_utilities.IsValueZero(surface_value, geometry_utilities.Tolerance1D()))
            marks[p] = SurfacePosition::OnSurface;
          else if (geometry_utilities.IsValueNegative(surface_value, geometry_utilities.Tolerance1D()))
            marks[p] = SurfacePosition::Negative;
          else
            throw std::runtime_error("unknown cell0D position");
        }

        return marks;
      };

      result.cell0Ds_surface_position = marking_function(result.cell0Ds_surface_value);
    }

    result.surface_intersections = IntersectMeshWithSurface(geometry_utilities,
                                                            result.cell0Ds_surface_value,
                                                            result.cell0Ds_surface_position,
                                                            mesh);

    result.surface_triangles = CreateTriangles(geometry_utilities,
                                               result.cell0Ds_surface_value,
                                               result.cell0Ds_surface_position,
                                               result.surface_intersections.surface_intersection_points,
                                               result.surface_intersections.cell0Ds_surface_points,
                                               result.surface_intersections.cell1Ds_surface_points,
                                               mesh);

    return result;
  }
  // ***************************************************************************
} // namespace Gedim
