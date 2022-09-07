/*!
 * @file graph.cpp
 *
 * @brief Graph class source file
 */

#include <assert.h>
#include <iostream>
#include <algorithm>
#include <stack>
#include <set>
#include <map>
#include <fstream>

#include "spoa/graph.hpp"

namespace spoa {

constexpr std::uint32_t kMaxAlphabetSize = 256;

std::unique_ptr<Node> Graph::createNode(std::uint32_t id, std::uint32_t code) {
    return std::unique_ptr<Node>(new Node(id, code));
}

Node::Node(std::uint32_t id, std::uint32_t code)
        : id_(id), code_(code), in_edges_(), out_edges_(),
        aligned_nodes_ids_() {
}

Node::~Node() {
}

bool Node::successor(std::uint32_t& dst, std::uint32_t label) const {

    for (const auto& edge: out_edges_) {
        for (const auto& l: edge->sequence_labels_) {
            if (l == label) {
                dst = edge->end_node_id_;
                return true;
            }
        }
    }
    return false;
}

std::shared_ptr<Edge> Node::get_edge(std::uint32_t end_id, bool is_outgoing_edge) const {
    if (is_outgoing_edge) {
        for (auto& edge: out_edges_) {
            if (edge->end_node_id_ == end_id) {
                return edge;
            }
        }
        return nullptr;
    } else {
        for (auto& edge: in_edges_) {
            if (edge->begin_node_id_ == end_id) {
                return edge;
            }
        }
        return nullptr;
    }
}

std::set<std::uint32_t> Node::get_associated_labels() const {
    std::set<std::uint32_t> label_set;
    for (const auto& edge: in_edges_) {
        for (const auto& label: edge->sequence_labels_) {
            label_set.insert(label);
        }
    }
    for (const auto& edge: out_edges_) {
        for (const auto& label: edge->sequence_labels_) {
            label_set.insert(label);
        }
    }
    return label_set;
}

std::unique_ptr<Edge> Graph::createEdge(std::uint32_t begin_node_id,
    std::uint32_t end_node_id, std::uint32_t label, std::uint32_t weight) {

    return std::unique_ptr<Edge>(new Edge(begin_node_id, end_node_id, label,
        weight));
}

Edge::Edge(std::uint32_t begin_node_id, std::uint32_t end_node_id,
    std::uint32_t label, std::uint32_t weight)
        : begin_node_id_(begin_node_id), end_node_id_(end_node_id),
        sequence_labels_(1, label), total_weight_(weight) {
}

Edge::~Edge() {
}

void Edge::add_sequence(std::uint32_t label, std::uint32_t weight) {
    sequence_labels_.emplace_back(label);
    total_weight_ += weight;
}

std::unique_ptr<Graph> createGraph() {
    return std::unique_ptr<Graph>(new Graph());
}

Graph::Graph()
        : num_sequences_(0), num_codes_(0), coder_(kMaxAlphabetSize, -1),
        decoder_(kMaxAlphabetSize, -1), nodes_(), rank_to_node_id_(),
        sequences_begin_nodes_ids_(), sequences_end_nodes_ids_(), consensus_() {
}

Graph::~Graph() {
}

std::uint32_t Graph::add_node(std::uint32_t code) {
    std::uint32_t node_id = nodes_.size();
    nodes_.emplace_back(createNode(node_id, code));
    return node_id;
}

void Graph::add_edge(std::uint32_t begin_node_id, std::uint32_t end_node_id,
    std::uint32_t weight) {

    assert(begin_node_id < nodes_.size() && end_node_id < nodes_.size());

    for (const auto& edge: nodes_[begin_node_id]->out_edges_) {
        if (edge->end_node_id_ == end_node_id) {
            edge->add_sequence(num_sequences_, weight);
            return;
        }
    }

    std::shared_ptr<Edge> edge = createEdge(begin_node_id, end_node_id,
        num_sequences_, weight);
    nodes_[begin_node_id]->out_edges_.emplace_back(edge);
    nodes_[end_node_id]->in_edges_.emplace_back(edge);
}

void Graph::add_edge(std::uint32_t begin_node_id, std::uint32_t end_node_id,
    std::uint32_t weight, std::uint32_t sequence_label) {

    assert(begin_node_id < nodes_.size() && end_node_id < nodes_.size());

    for (const auto& edge: nodes_[begin_node_id]->out_edges_) {
        if (edge->end_node_id_ == end_node_id) {
            edge->add_sequence(sequence_label, weight);
            return;
        }
    }

    std::shared_ptr<Edge> edge = createEdge(begin_node_id, end_node_id,
        sequence_label, weight);
    nodes_[begin_node_id]->out_edges_.emplace_back(edge);
    nodes_[end_node_id]->in_edges_.emplace_back(edge);
}

void Graph::add_alignment(const Alignment& alignment,
    const std::string& sequence, std::uint32_t weight) {

    add_alignment(alignment, sequence.c_str(), sequence.size(), weight);
}

void Graph::add_alignment(const Alignment& alignment, const char* sequence,
    std::uint32_t sequence_size, std::uint32_t weight) {

    std::vector<std::uint32_t> weights(sequence_size, weight);
    add_alignment(alignment, sequence, sequence_size, weights);
}

void Graph::add_alignment(const Alignment& alignment, const std::string& sequence,
    const std::string& quality) {

    add_alignment(alignment, sequence.c_str(), sequence.size(),
        quality.c_str(), quality.size());
}

void Graph::add_alignment(const Alignment& alignment, const char* sequence,
    std::uint32_t sequence_size, const char* quality,
    std::uint32_t quality_size) {

    std::vector<std::uint32_t> weights;
    for (std::uint32_t i = 0; i < quality_size; ++i) {
        weights.emplace_back(static_cast<std::uint32_t>(quality[i] - 32)); // PHRED quality
    }
    add_alignment(alignment, sequence, sequence_size, weights);
}

void Graph::add_alignment(const Alignment& alignment, const std::string& sequence,
    const std::vector<std::uint32_t>& weights) {

    add_alignment(alignment, sequence.c_str(), sequence.size(), weights);
}

void Graph::add_alignment(const Alignment& alignment, const char* sequence,
    std::uint32_t sequence_size, const std::vector<std::uint32_t>& weights) {

    if (sequence_size == 0) {
        ++num_sequences_;
        return;
    }
    if (sequence_size != weights.size()) {
        throw std::invalid_argument("[spoa::Graph::add_alignment] error: "
            "sequence and weights are of unequal size!");
    }

    for (std::uint32_t i = 0; i < sequence_size; ++i) {
        auto c = sequence[i];
        if (coder_[c] == -1) {
            coder_[c] = num_codes_;
            decoder_[num_codes_] = c;
            ++num_codes_;
        }
    }

    if (alignment.empty()) { // no alignment
        auto id_pairs = add_sequence(sequence, weights, 0, sequence_size);
        std::int32_t begin_node_id = id_pairs.first, end_node_id = id_pairs.second;
        ++num_sequences_;
        sequences_begin_nodes_ids_.emplace_back(begin_node_id);
        sequences_end_nodes_ids_.emplace_back(end_node_id);
        topological_sort();
        return;
    }

    std::vector<std::uint32_t> valid_seq_ids;
    for (const auto& it: alignment) {
        if (it.second != -1) {
            valid_seq_ids.emplace_back(it.second);
        }
    }

    assert(valid_seq_ids.front() <= sequence_size);
    assert(valid_seq_ids.back() + 1 <= sequence_size);

    std::uint32_t tmp = nodes_.size();
    std::int32_t begin_node_id = add_sequence(sequence, weights, 0,
        valid_seq_ids.front()).first;
    std::int32_t head_node_id = tmp == nodes_.size() ? -1 : nodes_.size() - 1;

    auto end_pair = add_sequence(sequence, weights, valid_seq_ids.back() + 1, sequence_size);
    std::int32_t tail_node_id = end_pair.first, end_node_id = end_pair.second;

    std::int32_t new_node_id = -1;
    float prev_weight = head_node_id == -1 ?
        0 : weights[valid_seq_ids.front() - 1];

    for (std::uint32_t i = 0; i < alignment.size(); ++i) {
        if (alignment[i].second == -1) {
            continue;
        }

        char letter = sequence[alignment[i].second];
        if (alignment[i].first == -1) {
            new_node_id = add_node(coder_[letter]);

        } else {
            if (decoder_[nodes_[alignment[i].first]->code_] == letter) {
                new_node_id = alignment[i].first;

            } else {
                std::int32_t aligned_to_node_id = -1;
                for (const auto& aid: nodes_[alignment[i].first]->aligned_nodes_ids_) {
                    if (decoder_[nodes_[aid]->code_] == letter) {
                        aligned_to_node_id = aid;
                        break;
                    }
                }

                if (aligned_to_node_id == -1) {
                    new_node_id = add_node(coder_[letter]);

                    for (const auto& aid: nodes_[alignment[i].first]->aligned_nodes_ids_) {
                        nodes_[new_node_id]->aligned_nodes_ids_.emplace_back(aid);
                        nodes_[aid]->aligned_nodes_ids_.emplace_back(new_node_id);
                    }

                    nodes_[new_node_id]->aligned_nodes_ids_.emplace_back(
                        alignment[i].first);
                    nodes_[alignment[i].first]->aligned_nodes_ids_.emplace_back(
                        new_node_id);

                } else {
                    new_node_id = aligned_to_node_id;
                }
            }
        }

        if (begin_node_id == -1) {
            begin_node_id = new_node_id;
        }

        if (head_node_id != -1) {
            // both nodes contribute to edge weight
            add_edge(head_node_id, new_node_id,
                prev_weight + weights[alignment[i].second]);
        }

        head_node_id = new_node_id;
        prev_weight = weights[alignment[i].second];
    }

    if (tail_node_id != -1) {
        // both nodes contribute to edge weight
        add_edge(head_node_id, tail_node_id,
            prev_weight + weights[valid_seq_ids.back() + 1]);
        sequences_end_nodes_ids_.emplace_back(end_node_id);
    } else {
        sequences_end_nodes_ids_.emplace_back(new_node_id);
    }

    ++num_sequences_;
    sequences_begin_nodes_ids_.emplace_back(begin_node_id);

    topological_sort();
}

std::pair<std::int32_t, std::int32_t> Graph::add_sequence(const char* sequence,
    const std::vector<std::uint32_t>& weights, std::uint32_t begin,
    std::uint32_t end) {

    if (begin == end) {
        return std::make_pair(-1, -1);
    }

    std::int32_t first_node_id = add_node(coder_[sequence[begin]]);

    std::uint32_t node_id;
    for (std::uint32_t i = begin + 1; i < end; ++i) {
        node_id = add_node(coder_[sequence[i]]);
        // both nodes contribute to edge weight
        add_edge(node_id - 1, node_id, weights[i - 1] + weights[i]);
    }

    return std::make_pair(first_node_id, node_id);
}

void Graph::topological_sort() {

    rank_to_node_id_.clear();

    // 0 - unmarked, 1 - temporarily marked, 2 - permanently marked
    std::vector<std::uint8_t> node_marks(nodes_.size(), 0);
    std::vector<bool> check_aligned_nodes(nodes_.size(), true);
    std::stack<std::uint32_t> nodes_to_visit;

    for (std::uint32_t i = 0; i < nodes_.size(); ++i) {
        if (node_marks[i] != 0) {
            continue;
        }

        nodes_to_visit.push(i);
        while (nodes_to_visit.size() != 0) {
            std::uint32_t node_id = nodes_to_visit.top();
            bool valid = true;

            if (node_marks[node_id] != 2) {
                for (const auto& edge: nodes_[node_id]->in_edges_) {
                    if (node_marks[edge->begin_node_id_] != 2) {
                        nodes_to_visit.push(edge->begin_node_id_);
                        valid = false;
                    }
                }

                if (check_aligned_nodes[node_id]) {
                    for (const auto& aid: nodes_[node_id]->aligned_nodes_ids_) {
                        if (node_marks[aid] != 2) {
                            nodes_to_visit.push(aid);
                            check_aligned_nodes[aid] = false;
                            valid = false;
                        }
                    }
                }

                assert((valid || node_marks[node_id] != 1) &&
                    "Graph is not a DAG!");

                if (valid) {
                    node_marks[node_id] = 2;
                    if (check_aligned_nodes[node_id]) {
                        rank_to_node_id_.push_back(node_id);
                        for (const auto& aid: nodes_[node_id]->aligned_nodes_ids_) {
                            rank_to_node_id_.emplace_back(aid);
                        }
                    }
                } else {
                    node_marks[node_id] = 1;
                }
            }

            if (valid) {
                nodes_to_visit.pop();
            }
        }
    }

    assert(is_topologically_sorted() == true);
}

bool Graph::is_topologically_sorted() const {
    assert(nodes_.size() == rank_to_node_id_.size());

    std::vector<bool> visited_nodes(nodes_.size(), false);
    for (std::uint32_t node_id: rank_to_node_id_) {
        for (const auto& edge: nodes_[node_id]->in_edges_) {
            if (visited_nodes[edge->begin_node_id_] == false) {
                return false;
            }
        }
        visited_nodes[node_id] = true;
    }

    return true;
}

std::uint32_t Graph::initialize_multiple_sequence_alignment(
    std::vector<std::uint32_t>& dst) const {

    dst.resize(nodes_.size(), 0);

    std::uint32_t msa_id = 0;
    for (std::uint32_t i = 0; i < nodes_.size(); ++i) {
        std::uint32_t node_id = rank_to_node_id_[i];

        dst[node_id] = msa_id;
        for (std::uint32_t j = 0; j < nodes_[node_id]->aligned_nodes_ids_.size(); ++j) {
            dst[rank_to_node_id_[++i]] = msa_id;
        }
        ++msa_id;
    }

    return msa_id;
}

void Graph::generate_multiple_sequence_alignment(std::vector<std::string>& dst,
    bool include_consensus) {

    // assign msa id to each node
    std::vector<std::uint32_t> node_id_to_msa_id;
    auto msa_length = initialize_multiple_sequence_alignment(node_id_to_msa_id);

    // extract sequences from graph and create msa strings (add indels(-) where
    // necessary)
    for (std::uint32_t i = 0; i < num_sequences_; ++i) {
        std::string alignment_str(msa_length, '-');
        std::uint32_t node_id = sequences_begin_nodes_ids_[i];

        while (true) {
            alignment_str[node_id_to_msa_id[node_id]] =
                decoder_[nodes_[node_id]->code_];

            if (!nodes_[node_id]->successor(node_id, i)) {
                break;
            }
        }

        dst.emplace_back(alignment_str);
    }

    if (include_consensus) {
        // do the same for consensus sequence
        traverse_heaviest_bundle_orig();

        std::string alignment_str(msa_length, '-');
        for (const auto& node_id: consensus_) {
            alignment_str[node_id_to_msa_id[node_id]] =
                decoder_[nodes_[node_id]->code_];
        }
        dst.emplace_back(alignment_str);
    }
}

void Graph::generate_multiple_sequence_alignment_custom(std::vector<std::string>& msa) {

    // assign msa id to each node
    std::vector<std::uint32_t> node_id_to_msa_id;
    auto msa_length = initialize_multiple_sequence_alignment(node_id_to_msa_id);

    // extract sequences from graph and create msa strings (add indels(-) where
    // necessary)
    // RITU: DO IT ONLY FOR DRAFT SEQ (SEQ NUMBER 0)
    //for (std::uint32_t i: interesting) {
    for (std::uint32_t i=0; i < 1; ++i) {
        std::string alignment_str(msa_length, '-');
        std::uint32_t node_id = sequences_begin_nodes_ids_[i];

        while (true) {
            alignment_str[node_id_to_msa_id[node_id]] =
                decoder_[nodes_[node_id]->code_];

            if (!nodes_[node_id]->successor(node_id, i)) {
                break;
            }
        }

        msa.emplace_back(alignment_str);
    }
    
    // do the same for consensus sequence
    traverse_heaviest_bundle_orig();

    std::string alignment_str(msa_length, '-');
    for (const auto& node_id: consensus_) {
        alignment_str[node_id_to_msa_id[node_id]] =
            decoder_[nodes_[node_id]->code_];
    }
    msa.emplace_back(alignment_str);
    
}

// Ritu
std::string Graph::generate_consensus_custom(std::vector<uint32_t>& dst) {

    auto consensus_str = this->generate_consensus();
   
    dst.resize(consensus_.size(), 0);

    std::vector<uint32_t> node_id_to_msa_id;
    initialize_multiple_sequence_alignment(node_id_to_msa_id);

    for (uint32_t i = 0; i < num_sequences_; ++i) {
        auto node_id = sequences_begin_nodes_ids_[i];

        uint32_t c = 0;
        while (true) {
            for (; c < consensus_.size() &&
                node_id_to_msa_id[consensus_[c]] < node_id_to_msa_id[node_id]; ++c);
            if (c >= consensus_.size()) {
                break;
            }

            if (node_id_to_msa_id[consensus_[c]] == node_id_to_msa_id[node_id]) {
                char letter = decoder_[nodes_[node_id]->code_];
                if (letter == consensus_str[c]) {
                    ++dst[c];
                }
            }

            if (!nodes_[node_id]->successor(node_id, i)) {
                break;
            }
        }
    }
    

    return consensus_str;
}

// Ritu
std::string Graph::generate_consensus_custom2(const std::vector<uint32_t>& interesting,std::vector<uint32_t>& dst) {

    auto consensus_str = this->generate_consensus();
   
    dst.resize(consensus_.size(), 0);

    std::vector<uint32_t> node_id_to_msa_id;
    initialize_multiple_sequence_alignment(node_id_to_msa_id);

    for (uint32_t i: interesting) {
        auto node_id = sequences_begin_nodes_ids_[i];

        uint32_t c = 0;
        while (true) {
            for (; c < consensus_.size() &&
                node_id_to_msa_id[consensus_[c]] < node_id_to_msa_id[node_id]; ++c);
            if (c >= consensus_.size()) {
                break;
            }

            if (node_id_to_msa_id[consensus_[c]] == node_id_to_msa_id[node_id]) {
                char letter = decoder_[nodes_[node_id]->code_];
                if (letter == consensus_str[c]) {
                    ++dst[c];
                }
            }

            if (!nodes_[node_id]->successor(node_id, i)) {
                break;
            }
        }
    }
    

    return consensus_str;
}

std::string Graph::generate_consensus() {
    exclude_labels_.clear();
    traverse_heaviest_bundle_orig();
    std::string consensus_str = "";
    for (const auto& node_id: consensus_) {
        consensus_str += decoder_[nodes_[node_id]->code_];
    }
    return consensus_str;
}

void Graph::generate_multiple_consensus(double threshold, std::vector<std::string>& all_consensuses) {
    exclude_labels_.clear();
    std::uint32_t C = 0;
    while (C < 1.0 * num_sequences_ * threshold) {
        C += traverse_heaviest_bundle();
        std::string consensus_str = "";
        for (const auto& node_id: consensus_) {
            consensus_str += decoder_[nodes_[node_id]->code_];
        }
        if (!all_consensuses.empty() && consensus_str == all_consensuses.back()) return;
        all_consensuses.push_back(consensus_str);
    }
}

void Graph::generate_multiple_consensus_new(double ratio_threshold, int constant_threshold, std::vector<std::string>& all_consensuses) {
    std::set<std::uint32_t> read_ids;
    for(std::uint32_t i = 0; i < num_sequences_; i++) read_ids.insert(i);
    std::unordered_set<std::string> result_consensus;
    recursive_consensus(ratio_threshold, constant_threshold, read_ids, result_consensus);
    for(auto & c : result_consensus) all_consensuses.push_back(c);
}

void Graph::recursive_consensus(const double ratio_threshold, const int constant_threshold, const std::set<std::uint32_t> & enabled_reads, std::unordered_set<std::string> & result_consensus) {
    std::set<std::uint32_t> bottleneck_reads;
    traverse_heaviest_bundle_new(enabled_reads, bottleneck_reads);
    
    if(bottleneck_reads.size() < constant_threshold || ((double)bottleneck_reads.size() / (double)enabled_reads.size()) > ratio_threshold) { // should reach this at some point
        std::string consensus_str = "";
        for (const auto& node_id: consensus_) {
            consensus_str += decoder_[nodes_[node_id]->code_];
        }
        result_consensus.insert(consensus_str);
    } else {        
        std::set<std::uint32_t> remainder;
        std::set_difference(enabled_reads.begin(), enabled_reads.end(), bottleneck_reads.begin(), bottleneck_reads.end(), std::inserter(remainder, remainder.end()));
        recursive_consensus(ratio_threshold, constant_threshold, bottleneck_reads, result_consensus);
        recursive_consensus(ratio_threshold, constant_threshold, remainder, result_consensus);
    }
}

std::pair<std::string, std::string> Graph::generate_consensus(double support_fraction) {
    exclude_labels_.clear();
    std::int32_t bottleneck = traverse_heaviest_bundle();
    std::string consensus_str = "";
    for (const auto& node_id: consensus_) {
        consensus_str += decoder_[nodes_[node_id]->code_];
    }
    
    if (bottleneck >= 1.0 * num_sequences_ * support_fraction) {
        return {consensus_str, ""};
    }
    traverse_heaviest_bundle();
    std::string consensus_str_2 = "";
    for (const auto& node_id: consensus_) {
        consensus_str_2 += decoder_[nodes_[node_id]->code_];
    }
    return {consensus_str, consensus_str_2};
}

std::string Graph::generate_consensus(std::vector<std::uint32_t>& dst,
    bool verbose) {

    auto consensus_str = generate_consensus();

    dst.clear();
    if (verbose == false) {
        for (const auto& node_id: consensus_) {
            std::uint32_t total_coverage = nodes_[node_id]->coverage();
            for (const auto& aid: nodes_[node_id]->aligned_nodes_ids_) {
                total_coverage += nodes_[aid]->coverage();
            }
            dst.emplace_back(total_coverage);
        }
    } else {
        dst.resize((num_codes_ + 1) * consensus_.size(), 0);

        std::vector<std::uint32_t> node_id_to_msa_id;
        initialize_multiple_sequence_alignment(node_id_to_msa_id);

        for (std::uint32_t i = 0; i < num_sequences_; ++i) {
            auto node_id = sequences_begin_nodes_ids_[i];

            bool count_indels = false;
            std::uint32_t c = 0, l;
            while (true) {
                for (; c < consensus_.size() &&
                    node_id_to_msa_id[consensus_[c]] < node_id_to_msa_id[node_id]; ++c);
                if (c >= consensus_.size()) {
                    break;
                }

                if (node_id_to_msa_id[consensus_[c]] == node_id_to_msa_id[node_id]) {
                    if (count_indels) {
                        for (std::uint32_t j = l + 1; j < c; ++j) {
                            ++dst[num_codes_ * consensus_.size() + j];
                        }
                    }
                    count_indels = true;
                    l = c;

                    ++dst[nodes_[node_id]->code_ * consensus_.size() + c];
                }

                if (!nodes_[node_id]->successor(node_id, i)) {
                    break;
                }
            }
        }
    }

    return consensus_str;
}



void find_mode(std::vector<std::uint32_t>& arr, std::set<std::uint32_t>& exclude, std::vector<std::uint32_t>& ans) {
    ans.clear();
    std::map<std::uint32_t, std::int32_t> counter;
    for (size_t i = 0; i < arr.size(); i++) {
        if (exclude.find(i) != exclude.end()) continue;
        auto x = arr[i];
        counter[x]++;
    }
    std::int32_t max_count = 0;
    for (const auto& pair: counter) {
        if (pair.second > max_count) {
            max_count = pair.second;
        }
    }
    for (const auto& pair: counter) {
        if (pair.second == max_count) {
            ans.push_back(pair.first);
        }
    }
}

void find_mode_new(std::vector<std::uint32_t>& arr, const std::set<std::uint32_t>& include, std::vector<std::uint32_t>& ans) {
    ans.clear();
    std::map<std::uint32_t, std::int32_t> counter;
    for (size_t i = 0; i < arr.size(); i++) {
        if (include.find(i) == include.end()) continue;
        auto x = arr[i];
        counter[x]++;
    }
    std::int32_t max_count = 0;
    for (const auto& pair: counter) {
        if (pair.second > max_count) {
            max_count = pair.second;
        }
    }
    for (const auto& pair: counter) {
        if (pair.second == max_count) {
            ans.push_back(pair.first);
        }
    }
}

void Graph::traverse_heaviest_bundle_orig() {

    std::vector<std::int32_t> predecessors(nodes_.size(), -1);
    std::vector<std::int32_t> predecessors_min_read_id(nodes_.size(), -1);
    std::vector<std::int64_t> scores(nodes_.size(), -1);

    std::uint32_t max_score_id = 0;
    std::int32_t current_minimum_id = 0;
    for (const auto& node_id: rank_to_node_id_) {
        for (const auto& edge: nodes_[node_id]->in_edges_) {
            if (scores[node_id] < edge->total_weight_ || (scores[node_id] == edge->total_weight_ && scores[predecessors[node_id]] < scores[edge->begin_node_id_])) {

                scores[node_id] = edge->total_weight_;
                predecessors[node_id] = edge->begin_node_id_;
                
                current_minimum_id = -1;
                for (auto lb: edge->sequence_labels_) {
                    if(current_minimum_id == -1 || current_minimum_id > lb) current_minimum_id = lb;
                }
                predecessors_min_read_id[node_id] = current_minimum_id;
                // std::cout << "Get minimum read id: " << current_minimum_id << std::endl;
            } else if(scores[node_id] == edge->total_weight_) { // compare minimum read ID?
                current_minimum_id = -1;
                for (auto lb: edge->sequence_labels_) {
                    if(current_minimum_id == -1 || current_minimum_id > lb) current_minimum_id = lb;
                }
                
                // std::cout << "Compare minimum read id: " << current_minimum_id << " " << predecessors_min_read_id[node_id] << std::endl;
                if(current_minimum_id < predecessors_min_read_id[node_id]) {
                    predecessors[node_id] = edge->begin_node_id_;
                    predecessors_min_read_id[node_id] = current_minimum_id;
                }
            }
        }

        if (predecessors[node_id] != -1) {
            scores[node_id] += scores[predecessors[node_id]];
        }

        if (scores[max_score_id] < scores[node_id]) {
            max_score_id = node_id;
        }
    }

    if (nodes_[max_score_id]->out_edges_.size() != 0) {

        std::vector<std::uint32_t> node_id_to_rank(nodes_.size(), 0);
        for (std::uint32_t i = 0; i < nodes_.size(); ++i) {
            node_id_to_rank[rank_to_node_id_[i]] = i;
        }

        while (nodes_[max_score_id]->out_edges_.size() != 0) {
            max_score_id = branch_completion(scores, predecessors,
                node_id_to_rank[max_score_id]);
        }
    }

    // traceback
    consensus_.clear();
    while (predecessors[max_score_id] != -1) {
        consensus_.emplace_back(max_score_id);
        max_score_id = predecessors[max_score_id];
    }
    consensus_.emplace_back(max_score_id);

    std::reverse(consensus_.begin(), consensus_.end());
}

std::int32_t Graph::traverse_heaviest_bundle() {
    std::vector<std::int32_t> predecessors(nodes_.size(), -1);
    std::vector<std::int32_t> scores(nodes_.size(), -1);

    for (const auto& node_id: rank_to_node_id_) {
        int bestScore = -1, bestNbScore = -1, pred_neighbour = -1;
        for (const auto& edge: nodes_[node_id]->in_edges_) {
            int weight = edge->sequence_labels_.size();
            int neighbourId = edge->begin_node_id_;
            if (!exclude_labels_.empty()) {
                weight = 0;
                for (auto lb: edge->sequence_labels_) {
                    if (exclude_labels_.find(lb) == exclude_labels_.end()) {
                        weight++;
                    }
                }
            }
            if (weight == 0) continue;
            if (weight > bestScore || (weight == bestScore && scores[neighbourId] > bestNbScore)) {
                bestScore = weight;
                bestNbScore = scores[neighbourId];
                pred_neighbour = neighbourId;
            }
        }
        if (bestScore == -1) {
            scores[node_id] = 0;
            scores[node_id] = -1;
        } else {
            scores[node_id] = bestScore + bestNbScore;
            predecessors[node_id] = pred_neighbour;
        }
    }
    std::int32_t max_score_id = -1;
    std::vector<std::uint32_t> modes;
    find_mode(sequences_end_nodes_ids_, exclude_labels_, modes);
    for (auto x: modes) {
        if (max_score_id == -1 || scores[x] > scores[max_score_id]) {
            max_score_id = x;
        }
    }
    find_mode(sequences_begin_nodes_ids_, exclude_labels_, modes);
    std::set<std::int32_t> end_ids;
    for (auto x: modes) {
        end_ids.insert(x);
    }

    // traceback
    std::int32_t cur_id = max_score_id, bottleneck = num_sequences_ + 1;
    std::shared_ptr<Edge> bottleneck_edge = nullptr;
    consensus_.clear();
    while (cur_id != -1 && end_ids.find(cur_id) == end_ids.end()) {
        consensus_.emplace_back(cur_id);
        if (predecessors[cur_id] != -1 && bottleneck > scores[cur_id] - scores[predecessors[cur_id]]) {
            bottleneck_edge = nodes_[cur_id]->get_edge(predecessors[cur_id], 0);
            bottleneck = scores[cur_id] - scores[predecessors[cur_id]];
        } else if (predecessors[cur_id] != -1 && bottleneck == scores[cur_id] - scores[predecessors[cur_id]]
                    && bottleneck_edge != nullptr 
                    && bottleneck_edge->sequence_labels_.size() < nodes_[cur_id]->get_edge(predecessors[cur_id], 0)->sequence_labels_.size()) {
            bottleneck_edge = nodes_[cur_id]->get_edge(predecessors[cur_id], 0);   
        }
        cur_id = predecessors[cur_id];
    }
    if (cur_id != -1) {
        consensus_.emplace_back(cur_id);
    }
    if (bottleneck_edge != nullptr) {
        for (auto& label: bottleneck_edge->sequence_labels_) {
            exclude_labels_.insert(label);
        }
    }
    
    std::reverse(consensus_.begin(), consensus_.end());
    return bottleneck;
}

std::int32_t Graph::traverse_heaviest_bundle_new(const std::set<std::uint32_t> & included_labels, std::set<std::uint32_t> & bottleneck_labels) {
    std::vector<std::int32_t> predecessors(nodes_.size(), -1);
    std::vector<std::int32_t> scores(nodes_.size(), -1);

    for (const auto& node_id: rank_to_node_id_) {
        int bestScore = -1, bestNbScore = -1, pred_neighbour = -1;
        for (const auto& edge: nodes_[node_id]->in_edges_) {
            int weight = edge->sequence_labels_.size();
            int neighbourId = edge->begin_node_id_;
            if (!included_labels.empty()) {
                weight = 0;
                for (auto lb: edge->sequence_labels_) {
                    if (included_labels.find(lb) != included_labels.end()) {
                        weight++;
                    }
                }
            }
            if (weight == 0) continue;
            if (weight > bestScore || (weight == bestScore && scores[neighbourId] > bestNbScore)) {
                bestScore = weight;
                bestNbScore = scores[neighbourId];
                pred_neighbour = neighbourId;
            }
        }
        if (bestScore == -1) {
            scores[node_id] = 0;
            scores[node_id] = -1;
        } else {
            scores[node_id] = bestScore + bestNbScore;
            predecessors[node_id] = pred_neighbour;
        }
    }
    std::int32_t max_score_id = -1;
    std::vector<std::uint32_t> modes;
    find_mode_new(sequences_end_nodes_ids_, included_labels, modes);
    for (auto x: modes) {
        if (max_score_id == -1 || scores[x] > scores[max_score_id]) {
            max_score_id = x;
        }
    }
    find_mode_new(sequences_begin_nodes_ids_, included_labels, modes);
    std::set<std::int32_t> end_ids;
    for (auto x: modes) {
        end_ids.insert(x);
    }

    // traceback
    std::int32_t cur_id = max_score_id, bottleneck = num_sequences_ + 1;
    std::shared_ptr<Edge> bottleneck_edge = nullptr;
    consensus_.clear();
    while (cur_id != -1 && end_ids.find(cur_id) == end_ids.end()) {
        consensus_.emplace_back(cur_id);
        if (predecessors[cur_id] != -1 && bottleneck > scores[cur_id] - scores[predecessors[cur_id]]) {
            bottleneck_edge = nodes_[cur_id]->get_edge(predecessors[cur_id], 0);
            bottleneck = scores[cur_id] - scores[predecessors[cur_id]];
        } else if (predecessors[cur_id] != -1 && bottleneck == scores[cur_id] - scores[predecessors[cur_id]]
                    && bottleneck_edge != nullptr 
                    && bottleneck_edge->sequence_labels_.size() < nodes_[cur_id]->get_edge(predecessors[cur_id], 0)->sequence_labels_.size()) {
            bottleneck_edge = nodes_[cur_id]->get_edge(predecessors[cur_id], 0);   
        }
        cur_id = predecessors[cur_id];
    }
    if (cur_id != -1) {
        consensus_.emplace_back(cur_id);
    }
    if (bottleneck_edge != nullptr) {
        for (auto& label: bottleneck_edge->sequence_labels_) {
            if(included_labels.find(label) != included_labels.end()) bottleneck_labels.insert(label);
        }
    }
    
    std::reverse(consensus_.begin(), consensus_.end());
    return bottleneck;
}

std::uint32_t Graph::branch_completion(std::vector<std::int64_t>& scores,
    std::vector<std::int32_t>& predecessors, std::uint32_t rank) {

    std::uint32_t node_id = rank_to_node_id_[rank];
    for (const auto& edge: nodes_[node_id]->out_edges_) {
        for (const auto& o_edge: nodes_[edge->end_node_id_]->in_edges_) {
            if (o_edge->begin_node_id_ != node_id) {
                scores[o_edge->begin_node_id_] = -1;
            }
        }
    }

    std::int64_t max_score = 0;
    std::uint32_t max_score_id = 0;
    for (std::uint32_t i = rank + 1; i < rank_to_node_id_.size(); ++i) {

        std::uint32_t node_id = rank_to_node_id_[i];
        scores[node_id] = -1;
        predecessors[node_id] = -1;

        for (const auto& edge: nodes_[node_id]->in_edges_) {
            if (scores[edge->begin_node_id_] == -1) {
                continue;
            }
            std::uint32_t weight = edge->sequence_labels_.size();

            if (scores[node_id] < weight ||
                (scores[node_id] == weight &&
                scores[predecessors[node_id]] <= scores[edge->begin_node_id_])) {

                scores[node_id] = weight;
                predecessors[node_id] = edge->begin_node_id_;
            }
        }

        if (predecessors[node_id] != -1) {
            scores[node_id] += scores[predecessors[node_id]];
        }

        if (max_score < scores[node_id]) {
            max_score = scores[node_id];
            max_score_id = node_id;
        }
    }

    return max_score_id;
}

// backtracing from right to left!
void Graph::extract_subgraph_nodes(std::vector<bool>& dst,
    std::uint32_t begin_node_id, std::uint32_t end_node_id) const {

    dst.resize(nodes_.size(), false);

    std::stack<std::uint32_t> nodes_to_visit;
    nodes_to_visit.push(begin_node_id);

    while (nodes_to_visit.size() != 0) {
        std::uint32_t node_id = nodes_to_visit.top();
        nodes_to_visit.pop();

        if (dst[node_id] == false && node_id >= end_node_id) {
            for (const auto& edge: nodes_[node_id]->in_edges_) {
                nodes_to_visit.push(edge->begin_node_id_);
            }
            for (const auto& aid: nodes_[node_id]->aligned_nodes_ids_) {
                nodes_to_visit.push(aid);
            }

            dst[node_id] = true;
        }
    }
}

std::unique_ptr<Graph> Graph::subgraph(std::uint32_t begin_node_id,
    std::uint32_t end_node_id,
    std::vector<std::int32_t>& subgraph_to_graph_mapping) const {

    std::vector<bool> is_subgraph_node;
    extract_subgraph_nodes(is_subgraph_node, end_node_id, begin_node_id);

    // init subgraph
    auto subgraph = std::unique_ptr<Graph>(new Graph());
    subgraph->num_sequences_ = num_sequences_;
    subgraph->num_codes_ = num_codes_;
    subgraph->coder_ = std::vector<std::int32_t>(coder_);
    subgraph->decoder_ = std::vector<std::int32_t>(decoder_);

    // create mapping from subgraph to graph and vice versa and add nodes to
    // subgraph
    subgraph_to_graph_mapping.resize(nodes_.size(), -1);
    std::vector<std::int32_t> graph_to_subgraph_mapping(nodes_.size(), -1);

    for (std::uint32_t i = 0; i < is_subgraph_node.size(); ++i) {
        if (is_subgraph_node[i] == false) {
            continue;
        }

        std::uint32_t subgraph_id = subgraph->add_node(nodes_[i]->code_);
        graph_to_subgraph_mapping[i] = subgraph_id;
        subgraph_to_graph_mapping[subgraph_id] = i;
    }

    // add edges and aligned nodes
    for (std::uint32_t i = 0; i < is_subgraph_node.size(); ++i) {
        if (is_subgraph_node[i] == false) {
            continue;
        }

        std::uint32_t subgraph_id = graph_to_subgraph_mapping[i];

        for (const auto& edge: nodes_[i]->in_edges_) {
            if (graph_to_subgraph_mapping[edge->begin_node_id_] == -1) {
                continue;
            }
            subgraph->add_edge(graph_to_subgraph_mapping[edge->begin_node_id_],
                subgraph_id, edge->total_weight_);
        }
        for (const auto& aid: nodes_[i]->aligned_nodes_ids_) {
            if (graph_to_subgraph_mapping[aid] == -1) {
                continue;
            }
            subgraph->nodes_[subgraph_id]->aligned_nodes_ids_.emplace_back(
                graph_to_subgraph_mapping[aid]);
        }
    }

    subgraph->topological_sort();

    return subgraph;
}

void Graph::update_alignment(Alignment& alignment,
    const std::vector<std::int32_t>& subgraph_to_graph_mapping) const {

    for (std::uint32_t i = 0; i < alignment.size(); ++i) {
        if (alignment[i].first != -1) {
            alignment[i].first = subgraph_to_graph_mapping[alignment[i].first];
        }
    }
}

void Graph::print_dot(const std::string& path) const {

    if (path.empty()) {
        return;
    }

    std::ofstream out(path);

    std::vector<std::int32_t> in_consensus(nodes_.size(), -1);
    std::int32_t rank = 0;
    for (const auto& id: consensus_) {
        in_consensus[id] = rank++;
    }

    out << "digraph " << num_sequences_ << " {" << std::endl;
    out << "    graph [rankdir=LR]" << std::endl;
    for (std::uint32_t i = 0; i < nodes_.size(); ++i) {
        out << "    " << i << " [label = \"" << i << " - ";
        out << static_cast<char>(decoder_[nodes_[i]->code_]) << "\"";
        if (in_consensus[i] != -1) {
            out << ", style=filled, fillcolor=goldenrod1";
        }
        out << "]" << std::endl;

        for (const auto& edge: nodes_[i]->out_edges_) {
            out << "    " << i << " -> " << edge->end_node_id_;
            out << " [label = \"" << edge->total_weight_ << "\"";
            if (in_consensus[i] + 1 == in_consensus[edge->end_node_id_]) {
                out << ", color=goldenrod1";
            }
            out << "]" << std::endl;
        }
        for (const auto& aid: nodes_[i]->aligned_nodes_ids_) {
            if (aid > i) {
                out << "    " << i << " -> " << aid;
                out << " [style = dotted, arrowhead = none]" << std::endl;
            }
        }
    }
    out << "}" << std::endl;

    out.close();
}

void Graph::convert_to_adjacent_list(std::vector<std::vector<std::pair<std::uint32_t, std::int32_t>>>& adj_list) {
    adj_list.assign(nodes_.size() + 2, std::vector<std::pair<std::uint32_t, std::int32_t>>());
    for (const auto node_id: rank_to_node_id_) {
        for (const auto neighbour: nodes_[node_id]->out_edges_) {
            std::uint32_t weight = neighbour->sequence_labels_.size();
            if (!exclude_labels_.empty()) {
                weight = 0;
                for (const auto label: neighbour->sequence_labels_) {
                    if (exclude_labels_.find(label) == exclude_labels_.end()) weight++;
                }
            }

            if (weight > 0) {
                adj_list[node_id].emplace_back(neighbour->end_node_id_, weight);
            }
        }
    }
    
    std::vector<int> count;
    count.assign(nodes_.size(), 0);
    for (size_t i = 0; i < sequences_begin_nodes_ids_.size(); i++) {
        if (exclude_labels_.empty() || exclude_labels_.find(i) == exclude_labels_.end()) {
            count[sequences_begin_nodes_ids_[i]]++;
        }
    }
    // Head id: nodes_.size(), tail id: nodes_.size() + 1
    for (size_t i = 0; i < count.size(); i++) {
        if (count[i] > 0) {
            adj_list[nodes_.size()].emplace_back(i, count[i]);
        }
    } 

    count.assign(nodes_.size(), 0);
    for (size_t i = 0; i < sequences_end_nodes_ids_.size(); i++) {
        if (exclude_labels_.empty() || exclude_labels_.find(i) == exclude_labels_.end()) {
            count[sequences_end_nodes_ids_[i]]++;
        }
    }
    // Head id: nodes_.size(), tail id: nodes_.size() + 1
    for (size_t i = 0; i < count.size(); i++) {
        if (count[i] > 0) {
            adj_list[i].emplace_back(nodes_.size() + 1, count[i]);
        }
    } 
}

std::string Graph::find_consensus() {
    int head_id = nodes_.size(), tail_id = nodes_.size() + 1;
    std::vector<std::vector<std::pair<std::uint32_t, std::int32_t>>> adj_list;
    convert_to_adjacent_list(adj_list);
    dijsktra(head_id, tail_id, adj_list);

    std::string consensus_str = "";
    for (const auto& node_id: consensus_) {
        consensus_str += decoder_[nodes_[node_id]->code_];
    }
    exclude_labels_.clear();
    return consensus_str;
}

std::pair<std::string, std::string> Graph::find_consensus(double support_fraction) {
    int head_id = nodes_.size(), tail_id = nodes_.size() + 1;
    std::vector<std::vector<std::pair<std::uint32_t, std::int32_t>>> adj_list;
    convert_to_adjacent_list(adj_list);
    std::int32_t bottleneck = dijsktra(head_id, tail_id, adj_list);

    // std::cout << "bottleneck " << bottleneck << std::endl;
    std::string consensus_str = "";
    for (const auto& node_id: consensus_) {
        if (node_id < nodes_.size())
        consensus_str += decoder_[nodes_[node_id]->code_];
    }
   
    if (bottleneck == -1 || bottleneck >= 1.0 * num_sequences_ * support_fraction) {
        return {consensus_str, ""};
    }

    convert_to_adjacent_list(adj_list);
    dijsktra(head_id, tail_id, adj_list);

    std::string consensus_str_2 = "";
    for (const auto& node_id: consensus_) {
        if (node_id < nodes_.size())
        consensus_str_2 += decoder_[nodes_[node_id]->code_];
    }
    exclude_labels_.clear();
    return {consensus_str, consensus_str_2};
}

inline std::int32_t get_edge_weight(std::uint32_t begin_id, std::uint32_t end_id, 
    std::vector<std::vector<std::pair<std::uint32_t, std::int32_t>>>& adj_list) {
    for (auto& nb: adj_list[begin_id]) {
        if (nb.first == end_id)
            return nb.second;
    }
    return 0;
}

std::int32_t Graph::dijsktra(std::uint32_t source, std::uint32_t destination,
                     std::vector<std::vector<std::pair<std::uint32_t, std::int32_t>>>& adj_list) {
    const std::int32_t INF = 10000000;
    const std::uint32_t N = adj_list.size();
    std::vector<std::int32_t> previous, width;
    std::vector<bool> used;
    previous.assign(N, -1);
    width.assign(N, -INF);
    used.assign(N, 0);
    width[source] = INF;
    std::priority_queue<std::pair<std::int32_t, std::uint32_t>> pq;
    pq.emplace(width[source], source);
    while (!pq.empty()) {
        std::uint32_t u = pq.top().second;
        pq.pop(); 
        used[u] = 0;
        for (auto neighbour: adj_list[u]) {
            std::uint32_t v = neighbour.first;
            std::int32_t weight = neighbour.second;
            if (width[v] > std::min(weight, width[u])) continue;
            std::int32_t alt = std::max(width[v], std::min(weight, width[u]));
            if (alt > width[v]) {
                width[v] = alt;
                previous[v] = u;
                if (!used[v]) {
                    used[v] = 1;
                    pq.emplace(width[v], v);
                }
            } else if (alt == width[v] && weight > get_edge_weight((std::uint32_t) previous[v], v, adj_list)) {
                previous[v] = u;
                if (!used[v]) {
                    used[v] = 1;
                    pq.emplace(width[v], v);
                }
            }
        }
    }

    consensus_.clear();
    if (width[destination] == -INF) {
        return -1;
    }
    std::int32_t bottleneck = width[destination];
    std::int32_t cur = destination;
    std::pair<std::int32_t, std::int32_t> bottleneck_edge_id = {-1, -1};
    while (cur != -1) {
        if (bottleneck_edge_id.first == -1 && previous[cur] != -1) {
            std::uint32_t begin = (std::uint32_t) previous[cur], end = (std::uint32_t) cur;
            if (get_edge_weight(begin, end, adj_list) == bottleneck) {
                bottleneck_edge_id.first = previous[cur];
                bottleneck_edge_id.second = cur;
            }
        }
        consensus_.push_back(cur);
        cur = previous[cur];
    }
    consensus_.pop_back();
    std::reverse(consensus_.begin(), consensus_.end());
    consensus_.pop_back();
    if (bottleneck_edge_id.first != -1 && bottleneck_edge_id.second != -1) {
        if ((std::uint32_t) bottleneck_edge_id.first == source) {
            for (auto lb: sequences_begin_nodes_ids_) exclude_labels_.insert(lb);
        } else if ((std::uint32_t) bottleneck_edge_id.second == destination) {
            for (auto lb: sequences_end_nodes_ids_) exclude_labels_.insert(lb);
        } else {
            auto bottleneck_edge = nodes_[bottleneck_edge_id.first]->get_edge(bottleneck_edge_id.second, 1);
            for (auto lb: bottleneck_edge->sequence_labels_) {
                exclude_labels_.insert(lb);
            }
        }
        return width[destination];
    }
    return -1;
}

void Graph::merge_homopolymers() {
    std::set<std::uint32_t> neighbors;
    for (const auto id: rank_to_node_id_) {
        if (nodes_[id]->out_edges().empty() || nodes_[id]->out_edges().back()->sequence_labels_.empty()) continue;
        neighbors.clear();
        for (auto edge: nodes_[id]->out_edges()) {
            neighbors.insert(edge->end_node_id());
        }
        for (auto edge: nodes_[id]->out_edges()) {
            auto nb_node_id = edge->end_node_id();
            if (nodes_[nb_node_id]->out_edges().empty() || nodes_[nb_node_id]->out_edges().back()->sequence_labels_.empty()) continue;
            for (auto nb_edge: nodes_[nb_node_id]->out_edges()) {
                auto nbnb_id = nb_edge->end_node_id();
                if (neighbors.find(nbnb_id) != neighbors.end() && nodes_[id]->code() != nodes_[nb_node_id]->code() && nodes_[nb_node_id]->code() == nodes_[nbnb_id]->code()) {
                    auto primary_edge = nodes_[id]->get_edge(nbnb_id, 1);
                    auto alter_edge = edge->sequence_labels_.size() < nb_edge->sequence_labels_.size() ? edge : nb_edge;
                    if (primary_edge->sequence_labels_.size() >= alter_edge->sequence_labels_.size()) {
                        for (auto label: alter_edge->sequence_labels_) {
                            primary_edge->add_sequence(label);
                        }
                        alter_edge->sequence_labels_.clear();
                    } else {
                        for (auto label: primary_edge->sequence_labels_) {
                            edge->add_sequence(label);
                            nb_edge->add_sequence(label);
                        }
                        primary_edge->sequence_labels_.clear();
                    }
                }
            }
        }
    }
}

void Graph::clear() {
    std::fill(coder_.begin(), coder_.end(), -1);
    std::fill(decoder_.begin(), decoder_.end(), -1);
    nodes_.clear();
    rank_to_node_id_.clear();
    sequences_begin_nodes_ids_.clear();
    sequences_end_nodes_ids_.clear();
    consensus_.clear();
}

}
