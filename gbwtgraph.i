%module gbwtgraph
%feature("flatnested", "1");
%feature("notabstract") GBWTGraph;

%include <std_pair.i>
%include <std_vector.i>
%include <std_string.i>
%include <stdint.i>

%inline {
 #include <typeinfo>

}


%inline{
typedef long long int nid_t;
}


%typemap(in) nid_t  {
        $1 = static_cast<uint64_t>(SWIG_convert_long ($input));
};


%typemap(out) nid_t  {
        $result = scheme_make_integer_value((int)$1);
};

%apply nid_t {gbwt::node_type, node_type, uint64_t}


%typemap(out)  std::vector<handle_t>
{

int vector_size=$1.size();
Scheme_Object* list_out[vector_size] ;

for(std::size_t i = 0; i < vector_size; ++i) {
        Scheme_Object* scheme_pair;
        uint64_t tmp = reinterpret_cast<uint64_t&>($1.at(i)) >> 1 ;
        bool tmp2 = reinterpret_cast<uint64_t&>($1.at(i)) & 1 ;
        Scheme_Object* orientation;
        Scheme_Object* scheme_integer = scheme_make_integer_value(tmp);

        if(tmp2)
        {
                orientation =  scheme_make_byte_string("-")  ;
        }
        else
        {

                orientation = scheme_make_byte_string("+") ;
        }
        scheme_pair = scheme_make_pair(orientation, scheme_integer) ;
        list_out[i] = scheme_pair;
}


Scheme_Object* scheme_result =  scheme_build_list ( vector_size , list_out);
$result = scheme_result;
};



%typemap(in) vector<long long int>  {
        Scheme_Object* the_list = $input;
        int lg = (int) scheme_list_length( $input);
        std::vector<long long int> *return_val ;
        long long int *tmp;
        for(int i = 0; i < lg ; ++i)
        {

                scheme_get_long_long_val(scheme_car(the_list), tmp  ) ;
                return_val->push_back (reinterpret_cast<long long int&>(*tmp));
                the_list = scheme_cdr(the_list) ;
        }

        $1 = *return_val ;
};







/* %typemap(out)  std::vector<handle_t> */
/* { */
/*  */
/* int vector_size=$1.size(); */
/* Scheme_Object* list_out[vector_size] ; */
/*  */
/* for(std::size_t i = 0; i < vector_size; ++i) { */
/*         Scheme_Object* scheme_pair; */
/*         uint64_t tmp = reinterpret_cast<uint64_t&>($1.at(i)) >> 1 ; */
/*         bool tmp2 = reinterpret_cast<uint64_t&>($1.at(i)) & 1 ; */
/*         Scheme_Object* orientation; */
/*         Scheme_Object* scheme_integer = scheme_make_integer_value(tmp); */
/*  */
/*         if(tmp2) */
/*         { */
/*                 orientation =  scheme_make_byte_string("-")  ; */
/*         } */
/*         else */
/*         { */
/*  */
/*                 orientation = scheme_make_byte_string("+") ; */
/*         } */
/*         scheme_pair = scheme_make_pair(orientation, scheme_integer) ; */
/*         list_out[i] = scheme_pair; */
/* } */
/*  */
/*  */
/* Scheme_Object* scheme_result =  scheme_build_list ( vector_size , list_out); */
/* $result = scheme_result; */
/* }; */










%inline %{
        typedef uint64_t node_type ;
        typedef uint64_t size_type ;
        typedef std::pair<size_type , size_type> range_type ;

%}  ;





%{
#include<handlegraph/types.hpp>
#include <handlegraph/util.hpp>
#include <gbwtgraph/gbwtgraph.h>
#include <gbwtgraph/gfa.h>
#include <gbwt/gbwt.h>
#include <gbwt/utils.h>
#include <stdio.h>
#include <stdint.h>
#include <iostream>
#include <sstream>
#include <random>

using namespace std;
using namespace gbwtgraph;
using namespace gbwt;

typedef long long int nid_t;

 
SearchState* convert_void (void *input ) {

return (SearchState*) input   ;


}                             ;


SearchState pointer_to_searchstate (void *input ) {
  SearchState *tmp   =   (SearchState*) input; 
  return (SearchState) *tmp ;

}                             




long long int test_vector (vector<nid_t> nodes)
{
return nodes[0];

};

vector<void*> vector_test (vector<void*> tmp) 
{
        return tmp;

} ;



typedef const handle_t& ref_handle ;

typedef gbwt::vector_type   vector_type ;


class GRAPH : public GBWTGraph{
        public:
                GRAPH(const GBWT& gbwt_index,  const SequenceSource& sequence_source):GBWTGraph::GBWTGraph( gbwt_index,  sequence_source)   {}

                bool follow_edges  (const handle_t& handle, bool go_left, const function<bool(const handle_t&)>& iteratee)
                {
                        return  this->follow_edges_impl (handle, go_left, iteratee)  ;
                }

                bool for_each_handle (const function<bool(const handle_t&)>& iteratee, bool parallel = false)
                {
                        return this->for_each_handle_impl(iteratee,  parallel) ;

                }


                SearchState extend(SearchState state, nid_t node) const {return this->index->extend(state, node);};


                gbwt::SearchState extend(SearchState state, handle_t  handle) const {return this->index->extend(state, handle_to_node(handle) );};


                void serialize_members(std::ostream& out) const{
                     this->serialize_members (out)  ;
                }

                // Underlying implementation to "deserialize" method.
                // Load the sequences from the istream.
                // Throws sdsl::simple_sds::InvalidData if sanity checks fail.
                // User must call set_gbwt() before using the graph.
                void deserialize_members(std::istream& in){
                this-> deserialize_members(in)      ;
                }


                std::vector<handle_t> collect_succesive_edges  (nid_t node){
                        std::vector<handle_t>* edges =  new (std::vector<handle_t>);
                        handle_t  handle = get_handle(node) ;
                        std::function<bool(const handle_t&)> iteratee =
                                [&edges](handle_t handle) {
                                        edges->push_back(handle) ;
                                        return true;
                                };
                        bool is_end = this->follow_edges(handle, false, iteratee);
                        return  *edges;
                }


                std::vector<handle_t> collect_succesive_edges (handle_t  handle){
                        std::vector<handle_t>* edges =  new (std::vector<handle_t>);
                        std::function<bool(const handle_t&)> iteratee =
                                [&edges](handle_t handle) {
                                        edges->push_back(handle) ;
                                        return true;
                                };
                        bool is_end = this->follow_edges(handle, false, iteratee);
                        return  *edges;
                }




                std::vector<gbwt::SearchState> collect_succesive_paths (handle_t  handle){
                         std::vector<gbwt::SearchState>* paths =  new (std::vector<gbwt::SearchState>);
                         gbwt::SearchState  initial_state  = this->get_state(handle)  ;
                         std::function<bool(const SearchState&)> iteratee =
                                [&paths](gbwt::SearchState state) {
                                        paths->push_back(state) ;
                                        return true;
                                };
                        bool is_end = this->follow_paths(initial_state ,  iteratee);
                        return  *paths;
                }




                std::vector<gbwt::SearchState> collect_succesive_paths (nid_t  node){
                        return this->collect_succesive_paths (this->node_to_handle(node));

                }

                std::vector<gbwt::SearchState> collect_succesive_paths (handle_t  handle , int number_of_paths ){
                        std::vector<gbwt::SearchState>* paths =  new (std::vector<gbwt::SearchState>);
                        if (number_of_paths < 1)
                                return *paths;
                        gbwt::SearchState  initial_state  = this->get_state(handle)  ;
                        int counter = number_of_paths ;
                        std::function<bool(const SearchState&)> iteratee =
                                [&paths,&counter](gbwt::SearchState state) {
                                        paths->push_back(state) ;
                                        if (counter==0)
                                                return false;
                                        counter--;
                                        return true;
                                };
                        bool is_end = this->follow_paths(initial_state ,  iteratee);
                        return  *paths;
                }


                std::vector<gbwt::SearchState> collect_succesive_paths (nid_t node, int number_of_paths ){
                        return this->collect_succesive_paths (this->node_to_handle(node), number_of_paths);
                }




               vector<void*> extend_to_valid_states (SearchState state){
                        long long int node_id  = state.node;
                        std::vector<void*>* edge_states =  new (vector<void*>);
                        handle_t handle = get_handle(node_id) ;
                        std::function<bool(const handle_t&)> iteratee =
                                [&edge_states,&state, *this](handle_t handle) {
                                        SearchState *tmp = new(SearchState);
                                        *tmp = state;
                                         SearchState  *tmp1 = new(SearchState);  
                                         *tmp1 = this->index->extend( *tmp , handle_to_node(handle)) ;    
                                         edge_states->push_back((void*)tmp1); 
                                        return true;
                                };
                        bool is_end = this->follow_edges(handle, false, iteratee);
                        return  *edge_states;
                }



               vector<void*> extend_to_valid_states (void* vstate){
                        SearchState state = * (SearchState *)vstate ;
                       long long int node_id  = state.node;
                        std::vector<void*>* edge_states =  new (vector<void*>);
                        handle_t handle = get_handle(node_id) ;
                        std::function<bool(const handle_t&)> iteratee =
                                [&edge_states,&state, *this](handle_t handle) {
                                        SearchState *tmp = new(SearchState);
                                        *tmp = state;
                                         SearchState  *tmp1 = new(SearchState);  
                                         *tmp1 = this->index->extend( *tmp , handle_to_node(handle)) ;    
                                         edge_states->push_back((void*)tmp1); 
                                        return true;
                                };
                        bool is_end = this->follow_edges(handle, false, iteratee);
                        return  *edge_states;
                }



};




GRAPH   gfa_to_gbwtgraph (char* file){
        std::string cppString = file ;
        std::pair<std::unique_ptr<gbwt::GBWT>, std::unique_ptr<gbwtgraph::SequenceSource>>  tmp = gbwtgraph::gfa_to_gbwt(cppString);
        GBWT* tmp_1 =  tmp.first.release();
        SequenceSource* tmp_2 =    tmp.second.release();
        GRAPH* the_graph = new GRAPH(*tmp_1, *tmp_2);
        return  *the_graph ;

}


GBWT gfa_to_gbwt (char* file)
{
        std::string cppString = file ;
        std::pair<std::unique_ptr<gbwt::GBWT>, std::unique_ptr<gbwtgraph::SequenceSource>>  tmp = gbwtgraph::gfa_to_gbwt(cppString);
        GBWT* tmp_1 =  tmp.first.release();
        return  *tmp_1  ;

}


vector<void *> sample_distribution (vector<float> weights, vector<void*> nodes, int number_of_draws){
        std::random_device rd;
        std::mt19937 gen(rd());
        std::discrete_distribution<> d(weights.begin(), weights.end());
        vector<void*> *sample_container = new (std::vector<void*>);
        for(int n=0; n < number_of_draws; ++n)
        {
                sample_container->push_back(nodes[d(gen)]);
        };

      return *sample_container;
}




%};


/* typedef std::pair<size_type , size_type> range_type ; */

%typemap (out) range_type // std::pair<nid_t, nid_t>
{
        Scheme_Object* start_node = scheme_make_integer_value($1.first);
        Scheme_Object* end_node = scheme_make_integer_value($1.second);
        Scheme_Object*  node_pair = scheme_make_pair(start_node, end_node) ;
        $result = node_pair;
};

%apply range_type {gbwt::range_type,  std::pair<size_type, size_type>  }

%inline {
typedef const handle_t& ref_handle ;
}

%typemap(out) handle_t
{
  Scheme_Object* scheme_pair;
  nid_t tmp = handlegraph::as_integer($1) >> 1 ;
  bool tmp2 = handlegraph::as_integer($1) & 1 ;
  Scheme_Object* orientation;
  Scheme_Object* scheme_integer = scheme_make_integer_value(tmp);
  if(tmp2)
  {
      orientation =  scheme_make_byte_string("-") ;}
  else
  {

  orientation = scheme_make_byte_string("+") ;}

        $result = scheme_make_pair(orientation, scheme_integer);
};

%inline{

        inline static handle_t pack(const uint64_t& number, const bool& bit) {
                assert(number < (0x1ULL << 63));
                return handlegraph::as_handle((number << 1) | (bit ? 1 : 0));
        }
}

%typemap(in) ref_handle
{



        handle_t tmp; //= handlegraph::as_handle((uint64_t) 4);
        uint64_t  node_id;
        Scheme_Object* obj2;
        stringstream orientation;
         if (SCHEME_PAIRP($input)){
                 Scheme_Object*  obj1 =  SCHEME_CAR ($input);
                 if ( SCHEME_STRINGP(obj1)  )
                 {
                         orientation << SCHEME_STR_VAL(obj1);
                 };

                 obj2 =   SCHEME_CDR ($input);
                 node_id = (uint64_t)SWIG_convert_long(obj2);


        };

                bool is_left = false  ;
                if (orientation.str() == "-")  {is_left = true;};
                tmp= pack(node_id, is_left);

        $1 =  &tmp;
};




%typemap (out)  std::pair<std::string, std::pair<nid_t, nid_t>>
{
        Scheme_Object* the_string = scheme_make_byte_string($1.first.c_str()) ;
        Scheme_Object* start_node = scheme_make_integer_value($1.second.first);
        Scheme_Object* end_node = scheme_make_integer_value($1.second.second);
        Scheme_Object*  node_pair = scheme_make_pair(start_node, end_node) ;
        Scheme_Object*  pair = scheme_make_pair(the_string, node_pair) ;
        $result = pair;
} ;


%typemap(in) vector<void*>  {
        Scheme_Object* the_list = $input;
        int lg = (int) scheme_list_length( $input);
        vector<void*> *return_val = new vector<void*>  ;
        void *tmp;

        for(int i = 0; i < lg ; ++i)
        {

                return_val->push_back ((void*) scheme_car(the_list));
                the_list = scheme_cdr(the_list) ;
        }
        $1 = *return_val ;
};


%typemap(out)  vector<void*>{
        int vector_size=$1.size();
        Scheme_Object* list_out[vector_size] ;

        int  the_type;
        for(std::size_t i = 0; i < vector_size; ++i) {
                void* scheme_list_value =  $1.at(i);
                the_type = (int) SCHEME_TYPE (  $1.at(i) );
                cout << "\n";
                cout << the_type ;
                if (( the_type == 47 )|| ( the_type == 55 ) || ( the_type == 62 ) || ( the_type == 59))
                {list_out[i] = scheme_list_value;}
                else
                {
                        list_out[i] = SWIG_NewPointerObj ( scheme_list_value, $descriptor(void*) , $owner); 
                } ;
        }
        Scheme_Object* scheme_result = scheme_build_list(vector_size , list_out);
        $result = scheme_result;
};


%typemap(out)  vector<SearchState*>{
        int vector_size=$1.size();
        Scheme_Object* list_out[vector_size] ;

        for(std::size_t i = 0; i < vector_size; ++i) {
                list_out[i] = SWIG_NewPointerObj ($1.at(i), $descriptor(SearchState) , $owner); 
        }
        Scheme_Object* scheme_result = scheme_build_list(vector_size , list_out);
        $result = scheme_result;
};






%rename(equals) operator==;
%rename(unequals) operator!=;
%include "/usr/local/include/handlegraph/types.hpp";





class GRAPH {
        public:
                GRAPH(const GBWT& gbwt_index,  const SequenceSource& sequence_source) ;
                ~GRAPH();

                /* std::vector<handle_t> get_all_handle_edges (void) ; */
                /*  */
                /* std::vector<handle_t> get_all_handles (void)  ; */


                bool follow_edges  (ref_handle handle, bool go_left, const function<bool(const handle_t&)>& iteratee);

                bool for_each_handle (const function<bool(const handle_t&)>& iteratee, bool parallel = false) ;
                /* handle_t get_handle(char *value); */


                bool has_node(nid_t node_id) const;



                // Look up the handle for the node with the given ID in the given orientation.
                handle_t get_handle(uint64_t  node_id, bool is_reverse = false) const;


                // Get the ID from a handle.
                nid_t get_id(ref_handle  handle) const;

                // Get the orientation of a handle.
                bool get_is_reverse(ref_handle  handle) const;

                // Invert the orientation of a handle (potentially without getting its ID).
                handle_t flip(ref_handle  handle) const;

                // Get the length of a node.
                size_t get_length(ref_handle  handle) const;

                // Get the sequence of a node, presented in the handle's local forward
                // orientation.
                std::string get_sequence(ref_handle  handle) const;

                // Returns one base of a handle's sequence, in the orientation of the
                // handle.
                char get_base(ref_handle  handle, size_t index) const;

                // Returns a substring of a handle's sequence, in the orientation of the
                // handle. If the indicated substring would extend beyond the end of the
                // handle's sequence, the return value is truncated to the sequence's end.
                std::string get_subsequence(ref_handle  handle, size_t index, size_t size) const;

                // Return the number of nodes in the graph.
                size_t get_node_count() const;

                nid_t min_node_id() const;

                nid_t max_node_id() const;





                size_t get_degree(ref_handle  handle, bool go_left) const;

                // Returns true if there is an edge that allows traversal from the left
                // handle to the right handle.
                bool has_edge(ref_handle  left, const handle_t& right) const;

                //------------------------------------------------------------------------------

                /*
                   SerializableHandleGraph interface.
                 */


                void set_gbwt(const gbwt::GBWT& gbwt_index);

                /* set_gbwt(const gbwt::GBWT& gbwt_index) */
                /* { */
                /*         this->index = &gbwt_index; */
                /*  */
                /*         if(!(this->index->bidirectional())) */
                /*         { */
                /*                 throw InvalidGBWT("GBWTGraph: The GBWT index must be bidirectional"); */
                /*         } */
                /* } */


                uint32_t get_magic_number() const;


                void serialize_members(std::ostream& out) const;

                void deserialize_members(std::istream& in);

                //------------------------------------------------------------------------------


                // Returns `true` if the graph contains a translation from node ids to segment names.
                bool has_segment_names() const;

                // Returns (GFA segment name, semiopen node id range) containing the handle.
                // If there is no such translation, returns ("id", (id, id + 1)).
                std::pair<std::string, std::pair<nid_t, nid_t>> get_segment(ref_handle  handle) const;

                // Returns (GFA segment name, starting offset in the same orientation) for the handle.
                // If there is no translation, returns ("id", 0).
                std::pair<std::string, size_t> get_segment_name_and_offset(ref_handle  handle) const;

                // Returns the name of the original GFA segment corresponding to the handle.
                // If there is no translation, returns the node id as a string.
                std::string get_segment_name(ref_handle  handle) const;

                // Returns the starting offset in the original GFA segment corresponding to the handle
                // in the same orientation as the handle.
                // If there is no translation, returns 0.
                size_t get_segment_offset(ref_handle  handle) const;

                // Calls `iteratee` with each segment name and the semiopen interval of node ids
                // corresponding to it. Stops early if the call returns `false`.
                // In GBWTGraph, the segments are visited in sorted order by node ids.
                void for_each_segment(const std::function<bool(const std::string&, std::pair<nid_t, nid_t>)>& iteratee) const;

                // Calls `iteratee` with each inter-segment edge and the corresponding segment names
                // in the canonical orientation. Stops early if the call returns `false`.
                void for_each_link(const std::function<bool(const edge_t&, const std::string&, const std::string&)>& iteratee) const;

                //------------------------------------------------------------------------------

                /*
                   GBWTGraph specific interface.
                 */


                // Serialize the the graph into the output stream in the simple-sds format.
                void simple_sds_serialize(std::ostream& out) const;

                // Deserialize or decompress the graph from the input stream and set the given
                // GBWT index. Note that the GBWT index is essential for loading the structure.
                // Throws sdsl::simple_sds::InvalidData if sanity checks fail and `InvalidGBWT`
                // if the GBWT index is not bidirectional.
                void simple_sds_load(std::istream& in, const gbwt::GBWT& gbwt_index);

                // Returns the size of the serialized structure in elements.
                size_t simple_sds_size() const;

                // Convert gbwt::node_type to handle_t.gbwt::
                static handle_t node_to_handle(gbwt::node_type node) { return handlegraph::as_handle(node); }

                // Convert handle_t to gbwt::node_type.
                static node_type handle_to_node(ref_handle  handle) { return handlegraph::as_integer(handle); }

                // Get node sequence as a pointer and length.
                view_type get_sequence_view(ref_handle  handle) const;

                // Determine if the node sequence starts with the given character.
                bool starts_with(ref_handle  handle, char c) const;

                // Determine if the node sequence ends with the given character.
                bool ends_with(ref_handle  handle, char c) const;

                // Convert handle_t to gbwt::SearchState.
                // Note that the state may be empty if the handle does not correspond to a real node.
                SearchState get_state(ref_handle  handle) const { return this->index->find(handle_to_node(handle)); }

                // Convert handle_t to gbwt::BidirectionalState.
                // Note that the state may be empty if the handle does not correspond to a real node.
                gbwt::BidirectionalState get_bd_state(ref_handle  handle) const { return this->index->bdFind(handle_to_node(handle)); }



                SearchState extend(SearchState state, nid_t node) const ;


                /* gbwt::SearchState extend(SearchState state, handle_t  handle) const; */

                // Get the search state corresponding to the vector of handles.
                gbwt::SearchState find(const std::vector<handle_t>& path) const;

                // Get the bidirectional search state corresponding to the vector of handles.
                gbwt::BidirectionalState bd_find(const std::vector<handle_t>& path) const;

                // Visit all successor states of this state and call iteratee for the state.
                // Stop and return false if the iteratee returns false.
                // Note that this does not visit empty successor states.
                bool follow_paths(gbwt::SearchState state, const std::function<bool(const gbwt::SearchState&)>& iteratee) const
                {
                        return this->follow_paths(this->get_single_cache(), state, iteratee);
                }

                // Visit all predecessor/successor states of this state and call iteratee for the state.
                // Stop and return false if the iteratee returns false.
                // Note that this does not visit empty predecessor/successor states.
                // Each state corresponds to a path. Going backward extends the path left, while going
                // extends it right. When going backward, the state is for the reverse orientation.
                bool follow_paths(gbwt::BidirectionalState state, bool backward,
                                const std::function<bool(const gbwt::BidirectionalState&)>& iteratee) const
                {
                        return this->follow_paths(this->get_single_cache(), state, backward, iteratee);
                }

                //------------------------------------------------------------------------------

                /*
                   Cached GBWTGraph specific interface. Each thread must use a separate cache.
                 */

                // Return a cache for the GBWT index. Note: The cache is not thread-safe.
                gbwt::CachedGBWT get_cache() const { return gbwt::CachedGBWT(*(this->index), false); }

                // Return a single-node cache for the GBWT index. Mostly for internal use.
                // Note: The cache is not thread-safe.
                gbwt::CachedGBWT get_single_cache() const { return gbwt::CachedGBWT(*(this->index), true); }

                // Convert handle_t to gbwt::SearchState.
                /* gbwt::SearchState get_state(const gbwt::CachedGBWT& cache, ref_handle  handle) const */
                /* { */
                /*         return cache.find(handle_to_node(handle)); */
                /* } */

                // Convert handle_t to gbwt::BidirectionalState.
                gbwt::BidirectionalState get_bd_state(const gbwt::CachedGBWT& cache, ref_handle  handle) const
                {
                        return cache.bdFind(handle_to_node(handle));
                }

                // Get the search state corresponding to the vector of handles.
                gbwt::SearchState find(const gbwt::CachedGBWT& cache, const std::vector<handle_t>& path) const;

                // Get the bidirectional search state corresponding to the vector of handles.
                gbwt::BidirectionalState bd_find(const gbwt::CachedGBWT& cache, const std::vector<handle_t>& path) const;





                // Visit all successor states of this state and call iteratee for the state.
                // Stop and return false if the iteratee returns false.
                // Note that the state may be empty if no path continues to that node.
                bool follow_paths(const gbwt::CachedGBWT& cache, gbwt::SearchState state,
                                const std::function<bool(const gbwt::SearchState&)>& iteratee) const;

                // Visit all predecessor/successor states of this state and call iteratee for the state.
                // Stop and return false if the iteratee returns false.
                // Note that the state may be empty if no path continues to that node.
                // Each state corresponds to a path. Going backward extends the path left, while going
                // extends it right. When going backward, the state is for the reverse orientation.
                bool follow_paths(const gbwt::CachedGBWT& cache, gbwt::BidirectionalState state, bool backward,
                                const std::function<bool(const gbwt::BidirectionalState&)>& iteratee) const;

                // Loop over all the handles to next/previous (right/left) nodes. Passes
                // them to a callback which returns false to stop iterating and true to
                // continue. Returns true if we finished and false if we stopped early.
                bool cached_follow_edges(const gbwt::CachedGBWT& cache, ref_handle  handle, bool go_left,
                                const std::function<bool(const handle_t&)>& iteratee) const;

                //------------------------------------------------------------------------------



                std::vector<handle_t> collect_succesive_edges  (nid_t node);


                std::vector<handle_t> collect_succesive_edges (handle_t  handle) ;


                std::vector<gbwt::SearchState> collect_succesive_paths (handle_t  handle);

                std::vector<gbwt::SearchState> collect_succesive_paths (nid_t  node) ;

                std::vector<gbwt::SearchState> collect_succesive_paths (handle_t  handle , int number_of_paths );


                std::vector<gbwt::SearchState> collect_succesive_paths (nid_t  node  , int number_of_paths );

                vector<void*> extend_to_valid_states (void* vstate);

               /* vector<SearchState*> extend_to_valid_states (void* vstate); */


};


GRAPH gfa_to_gbwtgraph (char* file);

GBWT gfa_to_gbwt (char* file) ;


SearchState* convert_void (void *input ) ;


SearchState pointer_to_searchstate (void *input ); 
//------------------------------------------------------------------------------

/*
   Traverse all haplotype-consistent windows in the graph and call lambda() for each window.
   Uses multiple threads, so the lambda should be thread-safe.
   A window starts with the sequence of a node and is followed by window_size - 1 bases
   from subsequent nodes. If no extensions are possible, a shorter substring of
   length >= window_size also qualifies as a window.
 */
void for_each_haplotype_window(const GBWTGraph& graph, size_t window_size,
                const std::function<void(const std::vector<handle_t>&, const std::string&)>& lambda,
                bool parallel);

//------------------------------------------------------------------------------


%template(LongPair) std::pair<std::uint64_t, std::uint64_t>;

struct SearchState
{
        node_type node;
        range_type range;
        SearchState() : node(ENDMARKER), range(Range::empty_range()) {}
        size_type size() const { return Range::length(this->range); }
        SearchState(node_type node_id, range_type offset_range) : node(node_id), range(offset_range) {}
        SearchState(node_type node_id, size_type sp, size_type ep) : node(node_id), range(sp, ep) {}
        bool empty() const { return Range::empty(this->range); }
        bool operator==(SearchState another) const { return (this->node == another.node && this->range == another.range); }
        bool operator!=(SearchState another) const { return (this->node != another.node || this->range != another.range); }

};




%inline {
        range_type test_range_type  (void) {
                range_type  the_range ;
                the_range.first = 3 ;
                the_range.second= 4 ;
                return the_range;
        };


        range_type SearchState_get_range (SearchState* state){
                std::pair<size_type, size_type> the_range =  state->range;
                range_type return_value;
                return_value.first = the_range.first;
                return_value.second= the_range.second;
                return return_value;
                }   ;





};


vector<void*> sample_distribution (vector<float> weights, vector<void*> nodes, int number_of_draws);

/* nid_t test_vector (vector<nid_t> nodes); */

/* vector<void*> vector_test (vector<void*> tmp) ; */


range_type SearchState_get_range (SearchState* state);

class GBWT {
public:

  GBWT();
  GBWT(const GBWT& source);

  /* explicit GBWT(const std::vector<GBWT>& sources); */

  void swap(GBWT& another);
  void resample(size_type sample_interval);

  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v = nullptr, std::string name = "") const;
  void load(std::istream& in);

  void simple_sds_serialize(std::ostream& out) const;
  void simple_sds_load(std::istream& in);
  size_t simple_sds_size() const;


  size_type size() const { return this->header.size; }
  bool empty() const { return (this->size() == 0); }
  size_type sequences() const { return this->header.sequences; }
  size_type sigma() const { return this->header.alphabet_size; }
  size_type effective() const { return this->header.alphabet_size - this->header.offset; }

  std::pair<size_type, size_type> runs() const;
  size_type samples() const { return this->da_samples.size(); }

  bool bidirectional() const { return this->header.get(GBWTHeader::FLAG_BIDIRECTIONAL); }

//------------------------------------------------------------------------------

  /*
    Metadata interface.
  */

  bool hasMetadata() const { return this->header.get(GBWTHeader::FLAG_METADATA); }
  void addMetadata() { this->header.set(GBWTHeader::FLAG_METADATA); }
  void clearMetadata() { this->metadata.clear(); this->header.unset(GBWTHeader::FLAG_METADATA); };

//------------------------------------------------------------------------------

  SearchState find(node_type node) const { return gbwt::find(*this, node); }
//all template classes are problematic

  /* template<class Iterator> */
  /* SearchState find(Iterator begin, Iterator end) const { return gbwt::find(*this, begin, end); } */

  SearchState prefix(node_type node) const { return gbwt::prefix(*this, node); }

  /* template<class Iterator> */
  /* SearchState prefix(Iterator begin, Iterator end) const { return gbwt::prefix(*this, begin, end); } */

  SearchState extend(SearchState state, node_type node) const { return gbwt::extend(*this, state, node); }

  /* template<class Iterator> */
  /* SearchState extend(SearchState state, Iterator begin, Iterator end) const { return gbwt::extend(*this, state, begin, end); } */
  /*  */
  size_type locate(node_type node, size_type i) const { return gbwt::locate(*this, edge_type(node, i)); }

  size_type locate(edge_type position) const { return gbwt::locate(*this, position); }

//problematic
  /* std::vector<size_type> locate(node_type node, range_type range) const { return this->locate(SearchState(node, range)); } */
  /* std::vector<size_type> locate(SearchState state) const; */
  /*  */
  /* vector_type extract(size_type sequence) const { return gbwt::extract(*this, sequence); } */
  vector_type extract(edge_type position) const { return gbwt::extract(*this, position); }
  /* vector_type extract(edge_type position, size_type max_length) const { return gbwt::extract(*this, position, max_length); } */

//endpr
//------------------------------------------------------------------------------

  /*
    Bidirectional search interface. The queries check that the parameters are valid.
    On error or failed search, the return value is an empty bidirectional search state.
  */

  BidirectionalState bdFind(node_type node) const { return gbwt::bdFind(*this, node); }

  BidirectionalState bdExtendForward(BidirectionalState state, node_type node) const { return gbwt::bdExtendForward(*this, state, node); }

  BidirectionalState bdExtendBackward(BidirectionalState state, node_type node) const { return gbwt::bdExtendBackward(*this, state, node); }

//------------------------------------------------------------------------------

  /*
    Low-level interface: Nodes. The interface assumes that node identifiers are valid,
    except in contains() / hasEdge(). This can be checked with contains().
  */

  bool contains(node_type node) const
  {
    return ((node < this->sigma() && node > this->header.offset) || node == ENDMARKER);
  }

  bool contains(edge_type position) const
  {
    return (this->contains(position.first) && position.second < this->nodeSize(position.first));
  }

  bool contains(SearchState state) const
  {
    return (this->contains(state.node) && !(state.empty()) && state.range.second < this->nodeSize(state.node));
  }

  bool hasEdge(node_type from, node_type to) const
  {
    return (this->contains(from) && this->record(from).hasEdge(to));
  }

  std::vector<edge_type> edges(node_type from) const
  {
    return this->record(from).outgoing;
  }

  node_type firstNode() const { return this->header.offset + 1; }
  comp_type toComp(node_type node) const { return (node == 0 ? node : node - this->header.offset); }
  node_type toNode(comp_type comp) const { return (comp == 0 ? comp : comp + this->header.offset); }

  size_type nodeSize(node_type node) const { return this->bwt.size(this->toComp(node)); }
  bool empty(node_type node) const { return this->bwt.empty(this->toComp(node)); }

//------------------------------------------------------------------------------

  /*
    Low-level interface: Navigation and searching. The interface assumes that node
    identifiers are valid. This can be checked with contains().
  */

  // On error: invalid_edge().
  edge_type LF(node_type from, size_type i) const
  {
    if(from == ENDMARKER) { return this->endmarker().LF(i); }
    return this->record(from).LF(i);
  }

  // On error: invalid_edge().
  edge_type LF(edge_type position) const
  {
    if(position.first == ENDMARKER) { return this->endmarker().LF(position.second); }
    return this->record(position.first).LF(position.second);
  }

  // On error: invalid_offset().
  size_type LF(node_type from, size_type i, node_type to) const
  {
    return this->record(from).LF(i, to);
  }

  // On error: invalid_offset().
  size_type LF(edge_type position, node_type to) const
  {
    return this->record(position.first).LF(position.second, to);
  }

  // On error: Range::empty_range().
  range_type LF(node_type from, range_type range, node_type to) const
  {
    return this->record(from).LF(range, to);
  }

  // On error: Range::empty_range().
  range_type LF(SearchState state, node_type to) const
  {
    return this->record(state.node).LF(state.range, to);
  }

  // On error: Range::empty_range().
  range_type bdLF(SearchState state, node_type to, size_type& reverse_offset) const
  {
    return this->record(state.node).bdLF(state.range, to, reverse_offset);
  }

//------------------------------------------------------------------------------

  /*
    Low-level interface: Sequences. The interface assumes that node identifiers are
    valid. This can be checked with contains().
  */

  // Starting position of the sequence or invalid_edge() if something fails.
  edge_type start(size_type sequence) const { return this->LF(ENDMARKER, sequence); }

  // Returns the sampled document identifier or invalid_sequence() if there is no sample.
  size_type tryLocate(node_type node, size_type i) const
  {
    return this->da_samples.tryLocate(this->toComp(node), i);
  }

  // Returns the sampled document identifier or invalid_sequence() if there is no sample.
  size_type tryLocate(edge_type position) const
  {
    return this->da_samples.tryLocate(this->toComp(position.first), position.second);
  }





} ;
