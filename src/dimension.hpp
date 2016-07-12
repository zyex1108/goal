#ifndef goal_dimension_hpp
#define goal_dimension_hpp

#include <Phalanx_DimTag.hpp>

struct Dim : public PHX::DimTag
{
  Dim() {}
  char const* name() const;
  static Dim const& tag();
};

struct Elem : public PHX::DimTag
{
  Elem() {}
  char const* name() const;
  static Elem const& tag();
};

struct Node : public PHX::DimTag
{
  Node() {}
  char const* name() const;
  static Node const& tag();
};

struct QP : public PHX::DimTag
{
  QP() {}
  char const* name() const;
  static QP const& tag();
};

struct Dummy : public PHX::DimTag
{
  Dummy() {}
  char const* name() const;
  static Dummy const& tag();
};

#endif
