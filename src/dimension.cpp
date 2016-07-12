#include "dimension.hpp"

char const* Dim::name() const
{static const char n[] = "Dim"; return n;}

Dim const& Dim::tag()
{static const Dim myself; return myself;}

char const* Elem::name() const
{static const char n[] = "Elem"; return n;}

Elem const& Elem::tag()
{static const Elem myself; return myself;}

char const* Node::name() const
{static const char n[] = "Node"; return n;}

Node const& Node::tag()
{static const Node myself; return myself;}

char const* QP::name() const
{static const char n[] = "QP"; return n;}

QP const& QP::tag()
{static const QP myself; return myself;}

char const* Dummy::name() const
{static const char n[] = "Dummy"; return n;}

Dummy const& Dummy::tag()
{static const Dummy myself; return myself;}
