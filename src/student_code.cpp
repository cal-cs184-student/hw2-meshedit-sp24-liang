#include "student_code.h"
#include "CGL/vector3D.h"
#include "halfEdgeMesh.h"
#include "mutablePriorityQueue.h"

using namespace std;

namespace CGL
{
  // search whether an element is in the vector
  bool search(std::vector<VertexIter> vec, VertexIter elementToFind){
    auto it = std::find(vec.begin(), vec.end(), elementToFind);

    // Check if the element was found
    if (it != vec.end()) {
        return true;
    } else {
        return false;
    }
  }

  /**
   * Evaluates one step of the de Casteljau's algorithm using the given points and
   * the scalar parameter t (class member).
   *
   * @param points A vector of points in 2D
   * @return A vector containing intermediate points or the final interpolated vector
   */
  std::vector<Vector2D> BezierCurve::evaluateStep(std::vector<Vector2D> const &points)
  { 
    // TODO Part 1.
    
    vector<Vector2D> intermediate_point;
    for (int i = 0; i < points.size()-1; i++){
      Vector2D point = (1 - t) * points[i] + t * points[i + 1];
      intermediate_point.push_back(point);
    }
    return intermediate_point;

    // return std::vector<Vector2D>();
  }

  /**
   * Evaluates one step of the de Casteljau's algorithm using the given points and
   * the scalar parameter t (function parameter).
   *
   * @param points    A vector of points in 3D
   * @param t         Scalar interpolation parameter
   * @return A vector containing intermediate points or the final interpolated vector
   */
  std::vector<Vector3D> BezierPatch::evaluateStep(std::vector<Vector3D> const &points, double t) const
  {
    // TODO Part 2.

    vector<Vector3D> intermediate_point;
    for (int i = 0; i < points.size()-1; i++){
      Vector3D point = (1 - t) * points[i] + t * points[i + 1];
      intermediate_point.push_back(point);
    }
    return intermediate_point;

    // return std::vector<Vector3D>();
  }

  /**
   * Fully evaluates de Casteljau's algorithm for a vector of points at scalar parameter t
   *
   * @param points    A vector of points in 3D
   * @param t         Scalar interpolation parameter
   * @return Final interpolated vector
   */
  Vector3D BezierPatch::evaluate1D(std::vector<Vector3D> const &points, double t) const
  {
    // TODO Part 2.

    vector<Vector3D> intermediate_point = points;
    while (intermediate_point.size() > 1){
      intermediate_point = evaluateStep(intermediate_point, t);
    }
    return intermediate_point[0];

    // return Vector3D();
  }

  /**
   * Evaluates the Bezier patch at parameter (u, v)
   *
   * @param u         Scalar interpolation parameter
   * @param v         Scalar interpolation parameter (along the other axis)
   * @return Final interpolated vector
   */
  Vector3D BezierPatch::evaluate(double u, double v) const 
  {  
    // TODO Part 2.

    Vector3D point, final_point; //final_point: final point on the surface
    vector<Vector3D> intermediate_point_u; //4 blue points on the grey curve
    for (int i = 0; i < controlPoints.size(); i++){
      point = evaluate1D(controlPoints[i],u);
      intermediate_point_u.push_back(point);
    }
    final_point = evaluate1D(intermediate_point_u,v);
    return final_point;

    // return Vector3D();
  }

  Vector3D Vertex::normal( void ) const
  {
    // TODO Part 3.
    // Returns an approximate unit normal at this vertex, computed by
    // taking the area-weighted average of the normals of neighboring
    // triangles, then normalizing.

    HalfedgeCIter h, init, h_twin;
    FaceCIter f0, f;
    Vector3D sum, weighted_normal;
    h = halfedge();      // get the outgoing half-edge of the vertex
    init = h;
    f0 = h->face();
    sum = f0->normal(); 
    do {
        h_twin = h->twin(); // get the opposite half-edge
        f = h_twin->face();
        weighted_normal = f->normal(); 
        sum += weighted_normal;
        h = h_twin->next();      // move to the next outgoing half-edge of the vertex
    } while(h != init);          // keep going until we are back where we were
    return sum.unit();    
    
    // return Vector3D();
  }

  EdgeIter HalfedgeMesh::flipEdge( EdgeIter e0 )
  {
    // TODO Part 4.
    // This method should flip the given edge and return an iterator to the flipped edge.

    if ( !e0->isBoundary() ){
      /****************  before ************/ 
      //halfedges
      HalfedgeIter h0 = e0->halfedge();
      HalfedgeIter h1 = h0->next();
      HalfedgeIter h2 = h1->next();
      HalfedgeIter h3 = h0->twin();
      HalfedgeIter h4 = h3->next();
      HalfedgeIter h5 = h4->next();
      HalfedgeIter h6 = h1->twin();
      HalfedgeIter h7 = h2->twin();
      HalfedgeIter h8 = h4->twin();
      HalfedgeIter h9 = h5->twin();
      //vertices
      VertexIter v0 = h0->vertex();
      VertexIter v1 = h3->vertex();
      VertexIter v2 = h2->vertex();
      VertexIter v3 = h5->vertex();
      //edges
      EdgeIter e1 = h1->edge();
      EdgeIter e2 = h2->edge();
      EdgeIter e3 = h4->edge();
      EdgeIter e4 = h5->edge();
      //faces
      FaceIter f0 = h0->face();
      FaceIter f1 = h3->face();

      /****************  after ************/ 
      //halfedges
      h0->setNeighbors(h1, h3, v3, e0, f0);
      h1->setNeighbors(h2, h7, v2, e2, f0);
      h2->setNeighbors(h0, h8, v0, e3, f0);
      h3->setNeighbors(h4, h0, v2, e0, f1);
      h4->setNeighbors(h5, h9, v3, e4, f1);
      h5->setNeighbors(h3, h6, v1, e1, f1);
      h6->setNeighbors(h6->next(), h5, v2, e1, h6->face());
      h7->setNeighbors(h7->next(), h1, v0, e2, h7->face());
      h8->setNeighbors(h8->next(), h2, v3, e3, h8->face());
      h9->setNeighbors(h9->next(), h4, v1, e4, h9->face());
      //vertices
      v0->halfedge() = h2;
      v1->halfedge() = h5;
      v2->halfedge() = h3;
      v3->halfedge() = h0;
      //edges
      e0->halfedge() = h0;
      e1->halfedge() = h5;
      e2->halfedge() = h1;
      e3->halfedge() = h2;
      e4->halfedge() = h4;
      //faces
      f0->halfedge() = h0;
      f1->halfedge() = h3;

      return e0;

    } else {
      return EdgeIter();
    }

  }

  VertexIter HalfedgeMesh::splitEdge( EdgeIter e0 )
  {
    // TODO Part 5.
    // This method should split the given edge and return an iterator to the newly inserted vertex.
    // The halfedge of this vertex should point along the edge that was split, rather than the new edges.

    if ( !e0->isBoundary() ){
      /****************  before ************/ 
      //halfedges
      HalfedgeIter h0 = e0->halfedge();
      HalfedgeIter h1 = h0->next();
      HalfedgeIter h2 = h1->next();
      HalfedgeIter h3 = h0->twin();
      HalfedgeIter h4 = h3->next();
      HalfedgeIter h5 = h4->next();
      HalfedgeIter h6 = h1->twin();
      HalfedgeIter h7 = h2->twin();
      HalfedgeIter h8 = h4->twin();
      HalfedgeIter h9 = h5->twin();
      //vertices
      VertexIter v0 = h0->vertex();
      VertexIter v1 = h3->vertex();
      VertexIter v2 = h2->vertex();
      VertexIter v3 = h5->vertex();
      //edges
      EdgeIter e1 = h1->edge();
      EdgeIter e2 = h2->edge();
      EdgeIter e3 = h4->edge();
      EdgeIter e4 = h5->edge();
      //faces
      FaceIter f0 = h0->face();
      FaceIter f1 = h3->face();

      /****************  after ************/ 
      ////create new halfedge, vertex, edge, face
      //new halfedge
      HalfedgeIter h10 = newHalfedge();
      h10->next() = h0;
      h10->twin() = h0;
      h10->vertex() = v0;
      h10->edge() = e0;
      h10->face() = f0;
      HalfedgeIter h11 = newHalfedge();
      h11->next() = h0;
      h11->twin() = h0;
      h11->vertex() = v0;
      h11->edge() = e0;
      h11->face() = f0;
      HalfedgeIter h12 = newHalfedge();
      h12->next() = h0;
      h12->twin() = h0;
      h12->vertex() = v0;
      h12->edge() = e0;
      h12->face() = f0;
      HalfedgeIter h13 = newHalfedge();
      h13->next() = h0;
      h13->twin() = h0;
      h13->vertex() = v0;
      h13->edge() = e0;
      h13->face() = f0;
      HalfedgeIter h14 = newHalfedge();
      h14->next() = h0;
      h14->twin() = h0;
      h14->vertex() = v0;
      h14->edge() = e0;
      h14->face() = f0;
      HalfedgeIter h15 = newHalfedge();
      h15->next() = h0;
      h15->twin() = h0;
      h15->vertex() = v0;
      h15->edge() = e0;
      h15->face() = f0;
      //new vertex
      VertexIter v4 = newVertex();
      v4->halfedge() = h0;
      v4->position = Vector3D();
      //new position for v4
      Vector3D p0 = v4->halfedge()->vertex()->position;
      Vector3D p1 = v4->halfedge()->twin()->vertex()->position;
      v4->position = 0.5*(p0 + p1);
      //new edge
      EdgeIter e5 = newEdge();
      e5->halfedge() = h0;
      EdgeIter e6 = newEdge();
      e6->halfedge() = h0;
      EdgeIter e7 = newEdge();
      e7->halfedge() = h0;
      //new face
      FaceIter f2 = newFace();
      f2->halfedge() = h0;
      FaceIter f3 = newFace();
      f3->halfedge() = h0;

      ////halfedges
      //setneighbors(next, twin, vertex, edge, face)
      h0->setNeighbors(h3, h15, v1, e0, f1);
      h1->setNeighbors(h14, h6, v1, e1, f0);
      h2->setNeighbors(h12, h7, v2, e2, f3);
      h3->setNeighbors(h5, h10, v4, e7, f1);
      h4->setNeighbors(h10, h8, v0, e3, f2);
      h5->setNeighbors(h0, h9, v3, e4, f1);
      h6->setNeighbors(h6->next(), h1, v2, e1, h6->face());
      h7->setNeighbors(h7->next(), h2, v0, e2, h7->face());
      h8->setNeighbors(h8->next(), h4, v3, e3, h8->face());
      h9->setNeighbors(h9->next(), h5, v1, e4, h9->face());
      h10->setNeighbors(h11, h3, v3, e7, f2);
      h11->setNeighbors(h4, h12, v4, e5, f2);
      h12->setNeighbors(h13, h11, v0, e5, f3);
      h13->setNeighbors(h2, h14, v4, e6, f3);
      h14->setNeighbors(h15, h13, v2, e6, f0);
      h15->setNeighbors(h1, h0, v4, e0, f0);
      ////vertices
      v0->halfedge() = h4;
      v1->halfedge() = h1;
      v2->halfedge() = h2;
      v3->halfedge() = h5;
      v4->halfedge() = h3;
      ////edges
      e0->halfedge() = h0;
      e1->halfedge() = h1;
      e2->halfedge() = h2;
      e3->halfedge() = h4;
      e4->halfedge() = h5;
      e5->halfedge() = h11;
      e6->halfedge() = h13;
      e7->halfedge() = h3;
      ////faces
      f0->halfedge() = h1;
      f1->halfedge() = h0;
      f2->halfedge() = h4;
      f3->halfedge() = h2;

      return v4;

    } else {
      return VertexIter();
    }
  }



  void MeshResampler::upsample( HalfedgeMesh& mesh )
  {

    // TODO Part 6.
    // This routine should increase the number of triangles in the mesh using Loop subdivision.
    // One possible solution is to break up the method as listed below
    // 1. Compute new positions for all the vertices in the input mesh, using the Loop subdivision rule,
    // and store them in Vertex::newPosition. At this point, we also want to mark each vertex as being
    // a vertex of the original mesh.
    
    // 2. Compute the updated vertex positions associated with edges, and store it in Edge::newPosition.
    
    // 3. Split every edge in the mesh, in any order. For future reference, we're also going to store some
    // information about which subdivide edges come from splitting an edge in the original mesh, and which edges
    // are new, by setting the flat Edge::isNew. Note that in this loop, we only want to iterate over edges of
    // the original mesh---otherwise, we'll end up splitting edges that we just split (and the loop will never end!)
    
    // 4. Flip any new edge that connects an old and new vertex.

    // 5. Copy the new vertex positions into final Vertex::position.

    // iterate over all edges in the mesh
    VertexIter new_vertex;
    vector<VertexIter> original_vertex_list, new_vertex_list;
    // //scan for the first time, store the original vertices
    // for (EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++) {
    //     if (!e->isBoundary()) {
    //         original_vertex_list.push_back(e->halfedge()->vertex());
    //     }
    // }

    // scan for the second time, do the split, store the newly-generated vertices by split
    for (EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++) {
        if (!e->isBoundary()) {
            new_vertex = mesh.splitEdge(e);
            new_vertex_list.push_back(new_vertex);
        }
    }

    //scan for the third time, do the flip, only flip edge that connects old vertex and new vertex
    bool flag_1, flag_2;
    EdgeIter flipped;
    for (EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++) {
        if (!e->isBoundary()) {
            VertexIter p0 = e->halfedge()->vertex();
            VertexIter p1 = e->halfedge()->twin()->vertex();
            flag_1 = search(new_vertex_list, p0);
            flag_2 = search(new_vertex_list, p1);
            if (flag_1 != flag_2){
              flipped = mesh.flipEdge(e);
            }
        }
    }

    //update position
    bool flag_3;
    Vector3D A, B, C, D;
    for (EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++) {
        //update one vertex, in every iteration
        VertexIter vertex = e->halfedge()->vertex();
        flag_3 = search(new_vertex_list, vertex);
        if (flag_3){      //can find it in the list，meaning the vertex is newly added
            if (!e->isBoundary()) {
              //calculate position for new vertex
              A = e->halfedge()->vertex()->position;
              B = e->halfedge()->twin()->vertex()->position;
              C = e->halfedge()->next()->next()->vertex()->position;
              D = e->halfedge()->twin()->next()->next()->vertex()->position;
              e->halfedge()->vertex()->position = 0.375 * (A+B) + 0.125*(C+D);
            } else{
              Vector3D A = e->halfedge()->vertex()->position;
              Vector3D B = e->halfedge()->next()->vertex()->position;
              e->halfedge()->vertex()->position = 0.5 * (A+B);
            }
        } else {  //can't find it in the list，meaning the vertex is old
            if (!e->isBoundary()){
              size_t n = vertex->degree();
              float u;
              if (n == 3){
                u = 0.1875;
              }else{
                u = 3/(8*n);
              }

              //find all neighbors, and then calculate sum of neighbor position
              Vector3D neighbor_position_sum(0., 0., 0. );
              HalfedgeCIter h_twin,h;
              Vector3D v;
              h = vertex->halfedge();      // get the outgoing half-edge of the vertex
              do {
                  h_twin = h->twin(); // get the opposite half-edge
                  v = h_twin->vertex()->position; 
                  neighbor_position_sum += v;     
                  h = h_twin->next();               // move to the next outgoing half-edge of the vertex
              } while(h != vertex->halfedge());          // keep going until we are back where we were
              e->halfedge()->vertex()->position = (1 - n * u) * (vertex->position) + u * neighbor_position_sum;
            }
            
        }
       
    }

    

  }
}
