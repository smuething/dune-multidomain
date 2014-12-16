#ifndef DUNE_MULTIDOMAIN_DOFMAPPER_HH
#define DUNE_MULTIDOMAIN_DOFMAPPER_HH

#include <algorithm>
#include <utility>
#include <vector>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/typeindex.hh>
#include <dune/grid/common/exceptions.hh>

namespace Dune {
namespace PDELab {
namespace MultiDomain {

template<typename GV>
class DOFMapper
{

public:
  typedef typename GV::template Codim<0>::EntityPointer ElementPointer;
  typedef typename GV::template Codim<0>::Entity Element;
  typedef typename GV::Intersection Intersection;
  typedef typename Element::Geometry::LocalCoordinate Coordinate;

  DOFMapper(const Intersection& is)
    : intersection_(is)
    , inside_(DOFStorageOnInside(is))
    , elementPointer_(inside_ ? is.inside() : is.outside())
  {}

  const Element& element() const
  {
    return *elementPointer_;
  }

  int mapSubIndex(int i, int codim) const
  {
    if (codim == 0)
      return inside_ ? intersection_.indexInInside() : intersection_.indexInOutside();
    if (codim == GV::Grid::dimension - 1)
      {
        Coordinate c = (inside_ ? intersection_.geometryInInside() : intersection_.geometryInOutside()).corner(i);
        const VertexList& vl = vertexList(element().type());
        typename VertexList::const_iterator it = std::find_if(vl.begin(),vl.end(),[c](Coordinate rc) -> typename GV::ctype { rc -= c; return rc.two_norm() < 1e-10;}); // TODO: don't use lambda expression!
        assert(it != vl.end());
        return it - vl.begin();
      }
    DUNE_THROW(GridError,"currently, only degrees of freedom on the intersection and the vertices are supported!");
  }

private:
  const Intersection& intersection_;
  const bool inside_;
  ElementPointer elementPointer_;

  typedef std::vector<Coordinate> VertexList;
  typedef std::vector<VertexList> VertexListContainer;

  static VertexListContainer _vertexListContainer;

  static bool DOFStorageOnInside(const Intersection& is)
  {
    if (is.conforming())
      return true;
    const int insideLevel = is.inside()->level();
    const int outsideLevel = is.outside()->level();
    assert(insideLevel != outsideLevel);
    return insideLevel > outsideLevel;
  }

  static const VertexList& vertexList(Dune::GeometryType gt)
  {
    const std::size_t i = LocalGeometryTypeIndex::index(gt);
    if (_vertexListContainer[i].empty())
      _vertexListContainer[i] = buildVertexList(gt);

    return _vertexListContainer[i];
  }

  static VertexList buildVertexList(Dune::GeometryType gt)
  {
    const Dune::ReferenceElement<typename GV::Grid::ctype,GV::Grid::dimension>& refEl =
      Dune::ReferenceElements<typename GV::Grid::ctype,GV::Grid::dimension>::general(gt);
    VertexList vl(refEl.size(GV::Grid::dimension));
    for (std::size_t i = 0; i < vl.size(); ++i)
      vl[i] = refEl.position(i,GV::Grid::dimension);
    return std::move(vl);
  }

};

template<typename GV>
typename DOFMapper<GV>::VertexListContainer DOFMapper<GV>::_vertexListContainer(LocalGeometryTypeIndex::size(GV::dimension));

} // namespace MultiDomain
} // namespace PDELab
} // namespace Dune

#endif // DUNE_MULTIDOMAIN_DOFMAPPER_HH
