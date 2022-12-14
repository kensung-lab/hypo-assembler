/*!
 * @file graph.hpp
 *
 * @brief Graph class header file
 */

#pragma once

#include <cstdint>
#include <memory>
#include <string>
#include <vector>
#include <utility>
#include <unordered_set>
#include <set>
#include <iostream>
#include <queue>

namespace spoa {

class Node;
class Edge;

class Graph;
std::unique_ptr<Graph> createGraph();

using Alignment = std::vector<std::pair<std::int32_t, std::int32_t>>;

class Graph {
public:
    ~Graph();
    std::set<std::uint32_t> labels_of_first_css_;

    const std::vector<std::unique_ptr<Node>>& nodes() const {
        return nodes_;
    }

    const std::vector<std::uint32_t>& rank_to_node_id() const {
        return rank_to_node_id_;
    }

    std::uint32_t num_codes() const {
        return num_codes_;
    };

    std::uint8_t coder(std::uint8_t c) const {
        return coder_[c];
    }

    std::uint8_t decoder(std::uint8_t code) const {
        return decoder_[code];
    }

    void add_alignment(const Alignment& alignment, const std::string& sequence,
        std::uint32_t weight = 1);

    void add_alignment(const Alignment& alignment, const char* sequence,
        std::uint32_t sequence_size, std::uint32_t weight = 1);

    void add_alignment(const Alignment& alignment, const std::string& sequence,
        const std::string& quality);

    void add_alignment(const Alignment& alignment, const char* sequence,
        std::uint32_t sequence_size, const char* quality,
        std::uint32_t quality_size);

    void add_alignment(const Alignment& alignment, const std::string& sequence,
        const std::vector<std::uint32_t>& weights);

    void add_alignment(const Alignment& alignment, const char* sequence,
        std::uint32_t sequence_size, const std::vector<std::uint32_t>& weights);

    void generate_multiple_sequence_alignment(std::vector<std::string>& dst,
        bool include_consensus = false);

    std::pair<std::string, std::string> generate_consensus(double support_fraction);

    std::string generate_consensus();

    std::pair<std::string, std::string> find_consensus(double support_fraction);

    std::string find_consensus();

    void merge_homopolymers();

    void generate_multiple_consensus(double threshold, std::vector<std::string>& all_consensuses);
    
    void generate_multiple_consensus_new(double ratio_threshold, int constant_threshold, std::vector<std::string>& all_consensuses);

    // returns  base coverages or complete summary matrix if verbose equals true
    std::string generate_consensus(std::vector<std::uint32_t>& dst,
        bool verbose = false);

    std::unique_ptr<Graph> subgraph(std::uint32_t begin_node_id,
        std::uint32_t end_node_id,
        std::vector<std::int32_t>& subgraph_to_graph_mapping) const;

    void update_alignment(Alignment& alignment,
        const std::vector<std::int32_t>& subgraph_to_graph_mapping) const;

    void print_dot(const std::string& path) const;

    void clear();

    friend std::unique_ptr<Graph> createGraph();

private:
    Graph();
    Graph(const Graph&) = delete;
    const Graph& operator=(const Graph&) = delete;

    static std::unique_ptr<Node> createNode(std::uint32_t id, std::uint32_t code);

    static std::unique_ptr<Edge> createEdge(std::uint32_t begin_node_id,
        std::uint32_t end_node_id, std::uint32_t label, std::uint32_t weight);

    std::uint32_t add_node(std::uint32_t code);

    void add_edge(std::uint32_t begin_node_id, std::uint32_t end_node_id,
        std::uint32_t weight);

    void add_edge(std::uint32_t begin_node_id, std::uint32_t end_node_id,
        std::uint32_t weight, std::uint32_t sequence_label);

    std::pair<std::int32_t, std::int32_t> add_sequence(const char* sequence,
        const std::vector<std::uint32_t>& weights, std::uint32_t begin,
        std::uint32_t end);

    void topological_sort();

    bool is_topologically_sorted() const;

    std::int32_t traverse_heaviest_bundle();
    void traverse_heaviest_bundle_orig();
    
    
    std::int32_t traverse_heaviest_bundle_new(const std::set<std::uint32_t> & included_labels, std::set<std::uint32_t> & bottleneck_labels);
    void recursive_consensus(const double ratio_threshold, const int constant_threshold, const std::set<std::uint32_t> & enabled_reads, std::unordered_set<std::string> & result_consensus);

    void convert_to_adjacent_list(std::vector<std::vector<std::pair<std::uint32_t, std::int32_t>>>& adj_list);

    std::int32_t dijsktra(std::uint32_t source, std::uint32_t destination,
        std::vector<std::vector<std::pair<std::uint32_t, std::int32_t>>>& adj_list);

    std::uint32_t branch_completion(std::vector<std::int64_t>& scores,
        std::vector<std::int32_t>& predecessors,
        std::uint32_t rank);

    void extract_subgraph_nodes(std::vector<bool>& dst,
        std::uint32_t current_node_id, std::uint32_t end_node_id) const;

    std::uint32_t initialize_multiple_sequence_alignment(
        std::vector<std::uint32_t>& dst) const;

    std::uint32_t num_sequences_;
    std::uint32_t num_codes_;
    std::vector<std::int32_t> coder_;
    std::vector<std::int32_t> decoder_;
    std::vector<std::unique_ptr<Node>> nodes_;
    std::vector<std::uint32_t> rank_to_node_id_;
    std::vector<std::uint32_t> sequences_begin_nodes_ids_;
    std::vector<std::uint32_t> sequences_end_nodes_ids_;
    std::vector<std::uint32_t> consensus_;
    std::set<std::uint32_t> exclude_labels_;
};

class Node {
public:
    ~Node();

    std::uint32_t id() const {
        return id_;
    }

    std::uint32_t code() const {
        return code_;
    }


    const std::vector<std::shared_ptr<Edge>>& in_edges() const {
        return in_edges_;
    }

    const std::vector<std::shared_ptr<Edge>>& out_edges() const {
        return out_edges_;
    }

    const std::vector<std::uint32_t>& aligned_nodes_ids() const {
        return aligned_nodes_ids_;
    }

    std::uint32_t coverage() const {
        return get_associated_labels().size();
    }

    bool successor(std::uint32_t& dst, std::uint32_t label) const;

    std::set<std::uint32_t> get_associated_labels() const;

    std::shared_ptr<Edge> get_edge(std::uint32_t end_id, bool is_outgoing_edge) const;

    friend Graph;

private:
    Node(std::uint32_t id, std::uint32_t code);
    Node(const Node&) = delete;
    const Node& operator=(const Node&) = delete;

    std::uint32_t id_;
    std::uint32_t code_;
    std::vector<std::shared_ptr<Edge>> in_edges_;
    std::vector<std::shared_ptr<Edge>> out_edges_;
    std::vector<std::uint32_t> aligned_nodes_ids_;
};

class Edge {
public:
    ~Edge();

    std::uint32_t begin_node_id() const {
        return begin_node_id_;
    }

    std::uint32_t end_node_id() const {
        return end_node_id_;
    }

    friend Graph;
    friend Node;

private:
    Edge(std::uint32_t begin_node_id, std::uint32_t end_node_id,
        std::uint32_t label, std::uint32_t weight);
    Edge(const Edge&) = delete;
    const Edge& operator=(const Edge&) = delete;

    void add_sequence(std::uint32_t label, std::uint32_t weight = 1);

    std::uint32_t begin_node_id_;
    std::uint32_t end_node_id_;
    std::vector<std::uint32_t> sequence_labels_;
    std::int64_t total_weight_;
};

}
