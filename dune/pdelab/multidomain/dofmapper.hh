#ifndef DUNE_MULTIDOMAIN_DOFMAPPER_HH
#define DUNE_MULTIDOMAIN_DOFMAPPER_HH

#include <vector>
#include <map>
#include <utility>

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
        typename VertexList::const_iterator it = std::find_if(vl.begin(),vl.end(),[c](Coordinate rc) { rc -= c; return rc.two_norm() < 1e-10;}); // TODO: don't use lambda expression!
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
  typedef std::map<Dune::GeometryType,VertexList> VertexListMap;

  static VertexListMap vertexListMap_;

  static bool DOFStorageOnInside(const Intersection& is)
  {
    if (is.conforming())
      return true;
    const int insideLevel = is.inside().level();
    const int outsideLevel = is.outside().level();
    assert(insideLevel != outsideLevel);
    return insideLevel > outsideLevel;
  }

  static const VertexList& vertexList(Dune::GeometryType gt)
  {
    typename VertexListMap::iterator it = vertexListMap_.find(gt);
    if (it != vertexListMap_.end())
      return it->second;
    vertexListMap_.insert({gt,buildVertexList(gt)});
    return vertexListMap_[gt];
  }

  static VertexList buildVertexList(Dune::GeometryType gt)
  {
    const Dune::GenericReferenceElement<typename GV::Grid::ctype,GV::Grid::dimension>& refEl =
      Dune::GenericReferenceElements<typename GV::Grid::ctype,GV::Grid::dimension>::general(gt);
    VertexList vl(refEl.size(GV::Grid::dimension));
    for (std::size_t i = 0; i < vl.size(); ++i)
      vl[i] = refEl.position(i,GV::Grid::dimension);
    return std::move(vl);
  }

};

template<typename GV>
typename DOFMapper<GV>::VertexListMap DOFMapper<GV>::vertexListMap_;

} // namespace MultiDomain
} // namespace PDELab
} // namespace Dune

#endif // DUNE_MULTIDOMAIN_DOFMAPPER_HH