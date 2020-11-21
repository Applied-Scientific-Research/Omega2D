/*
 * Random.cpp - Creates a random number generator fit to a specific range
 *
 * (c)2020 Applied Scientific Research, Inc.
 *         Mark J Stock <markjstock@gmail.com>
 *         Blake B Hillier <blakehillier@mac.com>
 */
#pragma once
#include <random>

struct RandomGenerator {
  static int instances;
  //static std::random_device m_rd; //Will be used to obtain a seed for the random number engine
  std::mt19937 m_gen; //Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<float> m_dist;
  float m_lb;
  float m_ub;
  RandomGenerator(float _lb, float _ub) : m_lb(_lb), m_ub(_ub) {
    instances++;
    m_gen = std::mt19937(instances);
    m_dist = std::uniform_real_distribution<float>(0.0, 1.0);
  }

  float operator()() { return m_lb + (m_ub-m_lb)*m_dist(m_gen); }
};

int RandomGenerator::instances = 0;
