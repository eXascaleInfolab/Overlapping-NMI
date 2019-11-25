/*
 *     Copyright: Aaron F. McDaid 2011 aaronmcdaid@gmail.com
 *     See our paper on arXiv: "Normalized Mutual Information to evaluate overlapping
 *     community finding algorithms" by Aaron F. McDaid, Derek Greene, Neil Hurley.
 *
 *     This program is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     This program is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#include <cmath>
#include <iostream>
#include <fstream>
#include <exception>
#include <limits>

#include <vector>
#include <string>
#include <set>
#include <map>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <getopt.h>

#include "cmdline.h"
#include "cnl_header_reader.hpp"
#include "aaron_utils.hpp"


using namespace std;

/*
 * Examples of args
 * -z
 *  -f 1,2,3
 *  -a (optional arg)
 */

//exceptions
struct MissingFile {};
struct EmptyFile {};

template <typename NodeId>
void onmi(const char * file1, const char * file2, const char * syncfile
, const bool do_omega_also, const bool allnmis, float membership=1.f);

static int global_verbose_flag = 0;

//! \brief Arguments parser
struct ArgParser: gengetopt_args_info {
	ArgParser(int argc, char **argv) {
		auto  err = cmdline_parser(argc, argv, this);
		if(err)
			throw std::invalid_argument("Arguments parsing failed: " + std::to_string(err));
	}

	~ArgParser() {
		cmdline_parser_free(this);
	}
};


int main(int argc, char **argv)
{
	ArgParser  args_info(argc, argv);

	if(args_info.inputs_num + args_info.sync_given != 2) {
		cmdline_parser_print_help();
		exit(1);
	}
	if(args_info.verbose_flag) {
		global_verbose_flag = 1;
	}
	const char *file1 = args_info.sync_given ? args_info.sync_arg : args_info.inputs[0];
	const char *file2 = args_info.inputs[1 - args_info.sync_given];
	if(args_info.textid_flag)
		onmi<string>(file1, file2, args_info.sync_arg, args_info.omega_flag
			, args_info.allnmis_flag, args_info.membership_arg);
	else onmi<uint32_t>(file1, file2, args_info.sync_arg, args_info.omega_flag
		, args_info.allnmis_flag, args_info.membership_arg);
}

//typedef string NodeId;

template <typename NodeId>
using Grouping = vector< set< NodeId > >;

template <typename NodeId>
struct NodeToGroup : public unordered_map< NodeId, set<int> >
{
	int sharedGroups(const NodeId n_, const NodeId m_) const {
		// PP2(n_,m_);
		static const set<int> emptySet;

		const set<int> * nGrps;
		if(this->count(n_)==1)
			nGrps = &(this->find(n_)->second);
		else
			nGrps = &emptySet;
		// PP(nGrps->size());

		const set<int> * mGrps;
		if(this->count(m_)==1)
			mGrps = &(this->find(m_)->second);
		else
			mGrps = &emptySet;
		// PP(mGrps->size());

		vector<int> intersection;
		set_intersection(
				nGrps->begin(), nGrps->end(),
				mGrps->begin(), mGrps->end(),
				back_inserter(intersection)
				);
		const int inter = intersection.size();
		// PP(inter);

		assert(inter <= (int)nGrps->size());
		assert(inter <= (int)mGrps->size());

		return inter;
	}
};

struct OverlapMatrix {
	map< pair<int,int> , int> om; // the pair is an ordered pair;  count
	int N;

	OverlapMatrix(): om(), N(0)  {}
};

template <typename NodeId>
NodeId to_id(string val)
{
	return stoul(val);
}

template <>
string to_id<string>(string val)
{
	// Skip share part if exists
	size_t iend = val.find(':');
	if(iend != string::npos)
		val.resize(iend);
	return val;
}

template <typename NodeId>
Grouping<NodeId> fileToSet(const char *fname, float membership=1.f
, unordered_set<NodeId> *nodes=nullptr) {
	Grouping<NodeId> ss;
	ifstream input(fname);
	unless(input.is_open())
		throw  MissingFile();

	string  line;
	size_t  clsnum = 0;  // The number of clusters
	size_t  ndsnum = 0;  // The number of nodes
	parseHeader(input, line, clsnum, ndsnum);
	const bool  estimNodes = !ndsnum;  // Is the number of nodes parsed or estimated
	const size_t  cmsbytes = ndsnum ? 0 : inputSize(input, fname);
	if(!ndsnum || !clsnum)
		estimateSizes(ndsnum, clsnum, cmsbytes, membership);

	if(clsnum)
		ss.reserve(clsnum);

	if(nodes && ndsnum) {
		nodes->clear();
		if(ndsnum)
			nodes->reserve(ndsnum);
	}

	// Note: the processing is started from the read line
	size_t  mbscnt = 0;  // The number of member nodes in all clusters, should be >= ndsnum
	do {
		// Skip empty lines and comments
		if(line.empty() || line[0] == '#')
			continue;

		typename Grouping<NodeId>::value_type cl;
		// Load cluster
		istringstream fields(line);
		string field;
		fields >> field;
		// Skip cluster id if specified
		if(field[field.length() - 1] == '>' && !(fields >> field))
			continue;
		do {
			// Skip commented tails
			if(field.empty() || field[0] == '#')
				break;
			// Note: this algorithm does not support node shares (will be skipped)
			cl.insert(to_id<NodeId>(field));
		} while(fields >> field);

		//forEach(const string &field, amd::rangeOverStream(fields, "\t ")) {
		if(!cl.empty()) {
			if(nodes)
				nodes->insert(cl.begin(), cl.end());
			mbscnt += cl.size();
			ss.push_back(move(cl));
		}
		//else cerr << "Warning: ignoring empty clusters in the file: " << fname << endl;
	} while(getline(input, line));

	// Rehash the nodes decreasing the allocated space and number of buckets
	// for the faster iterating if required
	if(nodes && nodes->size() < nodes->bucket_count() * nodes->max_load_factor() * 0.8f)
		nodes->reserve(nodes->size());
	ss.shrink_to_fit();  // Trim preallocated space for the clusters

	if(ndsnum && mbscnt != ndsnum) {
		if(!estimNodes && mbscnt < ndsnum)
			fprintf(stderr, "WARNING, Specified number of nodes is incorrect (specified: %lu"
				") in the header of '%s'. The actual number of members: %lu\n", ndsnum, fname, mbscnt);
		else printf("# Average%s membership in '%s': %.4G (%lu / %lu)\n"
			, estimNodes ? " estimated" : "", fname, float(mbscnt) / ndsnum, mbscnt, ndsnum);
	}

	return ss;
}

template <typename NodeId>
NodeToGroup<NodeId> nodeToGroup(const Grouping<NodeId> &g) {
	NodeToGroup<NodeId> n2g;
	for(int grpId = 0; grpId < (int)g.size(); grpId++) {
		const typename Grouping<NodeId>::value_type &grp = g.at(grpId);
		for(const NodeId &n : grp) {
			n2g[n].insert(grpId);
		}
	}
	return n2g;
}

template <typename NodeId>
const OverlapMatrix overlapMatrix(const NodeToGroup<NodeId> &ng1, const NodeToGroup<NodeId> &ng2) {
	OverlapMatrix om;
	unordered_set< NodeId > nodes;
	for(const typename NodeToGroup<NodeId>::value_type &n : ng2)
		nodes.insert(n.first);
	for(const typename NodeToGroup<NodeId>::value_type &n: ng1)
		nodes.insert(n.first);
	for(auto const &n: nodes) {
		if(ng1.count(n)) for(const int g1: ng1.find(n)->second) {
			if(ng2.count(n)) for(const int g2: ng2.find(n)->second) {
				om.om[make_pair(g1,g2)] ++;
			}
		}
	}
	om.N=nodes.size();
	return om;
}

double H(const int x, const int N) {
	if(x==0)
		return 0.0;
	const double Px = double(x) / double(N);
	assert(x>0 && Px > 0.0);
	return x != N ? -x*log2(Px) : 1;  // Matching of the same clusters should yield 1
}

double h (const double p) {
	if(p < numeric_limits<double>::epsilon()) {
		assert(p >= 0);
		return 0;
	}
	return p <= 1 - numeric_limits<double>::epsilon() ? -p*log2(p) : 1;  // Matching of the same clusters should yield 1
}

double H_X_given_Y (const int y, const int x, const int o, const int N) {
	// the NON-NORMALIZED mutual information
	// given sets of size 'l' and 'r', where there are 'o' nodes in common, what's their similarity score?
		assert(o>0 && y>=o && x>=o && y <= N && x <= N);
		const double H_Y = H(y,N) + H(N-y,N);
		const double H_X = H(x,N) + H(N-x,N); // just used in the assertion

		const double Px0y0 = (N-x-y+o)  /double(N); // TODO delete these, and following, lines
		const double Px1y0 = (x-o)      /double(N);
		const double Px0y1 = (y-o)      /double(N);
		const double Px1y1 =  o         /double(N);
		DYINGWORDS(VERYCLOSE(Px0y0 + Px0y1 + Px1y0 + Px1y1, 1.0)) {
			PP(Px0y0);
			PP(Px1y0);
			PP(Px0y1);
			PP(Px1y1);
			PP(Px0y0 + Px0y1 + Px1y0 + Px1y1 - 1.0);
		}
		// PP(h(Px0y0));
		// PP(h(Px1y1));
		// PP(h(Px0y1));
		// PP(h(Px1y0));
		// PP( h(Px0y0)+h(Px1y1) -h(Px0y1)-h(Px0y1) );
		if( h(Px0y0)+h(Px1y1) <= h(Px0y1)+h(Px1y0) )
			return H_X;
		const double H_XY =
			  H(N-x-y+o, N)
			+ H(x-o, N)
			+ H(y-o, N)
			+ H(o, N)
			;
		if(H_X < H_XY - H_Y) { // this shouldn't really happen, except by a roundoff error
			DYINGWORDS( VERYCLOSE(H_X , H_XY - H_Y)) {
				PP(H_X);
				PP(H_XY - H_Y);
				PP(H_X - (H_XY - H_Y));
			}
			return H_X;
		} else
			return H_XY - H_Y;
}

template<typename NodeId>
double HX_given_BestY (const OverlapMatrix &om, const Grouping<NodeId> &g1, const Grouping<NodeId> &g2, const int realxId) {
	// realxId is a community in g2.
	// X is g2, Y is g1
	// we're looking for the bits to encode X_realxId given all of Y
	const int sizeOfXComm = g2.at(realxId).size();
	double bestSoFar = H(sizeOfXComm,om.N) + H(om.N-sizeOfXComm,om.N);
	map< pair<int,int> ,int >::const_iterator            i = om.om.lower_bound(make_pair(realxId  , numeric_limits<int>::min()));
	map< pair<int,int> ,int >::const_iterator   endOfRange = om.om.lower_bound(make_pair(realxId+1, numeric_limits<int>::min()));
	for(; i != endOfRange; ++           i)
	{
		int xId =  i->first.first;
		int yId =  i->first.second;
		assert(realxId == xId);
		const int overlap = i->second;
		const double H_XgivenY = H_X_given_Y(g1.at(yId).size(),g2.at(xId).size(), overlap, om.N);
		if(bestSoFar > H_XgivenY)
			bestSoFar = H_XgivenY;
		// cout << '\t'; PPt(g1.at(fromId).size());
		// PPt(g2.at(xId).size());
		// PPt(overlap);
		// PP(H_XgivenY);
	}
	return bestSoFar;
}

template<bool normalizeTooSoon, typename NodeId>
double VI_oneSide (const OverlapMatrix &om, const Grouping<NodeId> &g1, const Grouping<NodeId> &g2) {
	// this doesn't return the (N)MI. It's the non-mutual information, optionally normalized too soon
	const int N = om.N;
	double total = 0.0;
	for(int toId = 0; toId < (int)g2.size(); toId++) {
		const double unnorm = HX_given_BestY(om, g1, g2, toId);
		if(normalizeTooSoon) {
			const int x = g2.at(toId).size();
			const double H_X = H(x,N) + H(N-x,N);
			const double norm = unnorm / H_X; // might be NaN
			if(H_X < numeric_limits<double>::epsilon()) { // the communities take up the whole set of nodes, and hence won't need any bits to be encoded. No need to add anything to 'total'
				assert(unnorm < numeric_limits<double>::epsilon()
					&& "unnorm should be ~=0 when H_X ~= 0");

				// in this case norm is 0/0, but we'll just define this as 1 // This is the bugfix/ambiguityfix to make it the same as the LFK software
				total += 1.0;
			} else {
				unless(norm <= 1.01) {
					PP(x);
					PP(N);
					PP(H_X);
					PP(unnorm);
					PP(unnorm / H_X);
					PP(unnorm - H_X);
				}
				assert(norm <= 1.01 && norm >= -0.01
					&& "norm should E [0, 1]");
				total += norm;
			}
		} else {
			total += unnorm;
		}
	}
	if(normalizeTooSoon) {
		return total / g2.size(); // why is this total zero? When it should be one in some cases.
	} else
		return total;
}

template<typename NodeId>
double LFKNMI(const OverlapMatrix &om, const OverlapMatrix &omFlipped
, const Grouping<NodeId> &g1, const Grouping<NodeId> &g2) {
	return 1.0 - 0.5 *
		( VI_oneSide<true>(omFlipped, g1, g2)
		+ VI_oneSide<true>(om       , g2, g1) );
}

struct Max {
	double operator() (const double H_Xs, const double H_Ys) const {
		return H_Xs > H_Ys ? H_Xs : H_Ys;
	}
};

struct Sqrt {
	double operator() (const double H_Xs, const double H_Ys) const {
		return sqrt(H_Xs * H_Ys);
	}
};

struct Sum {
	double operator() (const double H_Xs, const double H_Ys) const {
		return 0.5 * (H_Xs + H_Ys);
	}
};

struct Min {
	double operator() (const double H_Xs, const double H_Ys) const {
		return H_Xs > H_Ys ? H_Ys : H_Xs;
	}
};

template<typename Combiner, typename NodeId>
double aaronNMI(const OverlapMatrix &om, const OverlapMatrix &omFlipped
, const Grouping<NodeId> &g1, const Grouping<NodeId> &g2) {
	double H_Xs = 0.0;
	for(int toId = 0; toId < (int)g2.size(); toId++) {
		const int x = g2.at(toId).size();
		H_Xs += H(x, om.N)+H(om.N-x, om.N);
	}
	double H_Ys = 0.0;
	for(int fromId = 0; fromId < (int)g1.size(); fromId++) {
		const int x = g1.at(fromId).size();
		H_Ys += H(x, om.N)+H(om.N-x, om.N);
	}
#ifdef DEBUG
	printf("aaronNMI(), H_Xs: %G (g2 size: %lu), H_Ys: %G (g1 size: %lu), VI1: %G, VI2: %G\n"
			, H_Xs, g2.size(), H_Ys, g1.size(), VI_oneSide<false>(omFlipped, g1, g2), VI_oneSide<false>(om, g2, g1));
#endif // DEBUG
	//return 0.5*( H_Xs+H_Ys - VI_oneSide<false>(omFlipped, g1, g2)
	//- VI_oneSide<false>(om, g2, g1) ) / Combiner()(H_Xs, H_Ys);
	// Note: initial formula does not work for the sqrt and nax combiners for the fully overlapping clusters:
	// the ground-truth can be types / multiple categories for the same items (cluster: 1 2; gt: 1 2; 1 2)
	// Internal Note: the counter example for this measure is rotshift_single.cnl
	return Combiner()(H_Xs - VI_oneSide<false>(omFlipped, g1, g2),
		H_Ys - VI_oneSide<false>(om, g2, g1)) / Combiner()(H_Xs, H_Ys);
}

typedef long long int lli;

template <typename NodeId>
pair<double,double> omega(const NodeToGroup<NodeId> &ng1, const NodeToGroup<NodeId> &ng2) {
	// To understand this implement this implementation, look at the formulation
	// on page 6 of this: http://www.aclweb.org/anthology/W12-2602
	// That defines:
	//
	//                      Obs(s1,s2) - Exp(s1,s2)
	//      Omega(s1,s2) =  -----------------------
	//                               1 - Exp(s1,s2)
	//
	// The only 'trick' I have used is to multiply above and below the line by N^2
	//
	// Finally, note that the Omega index can be negative.
	// Obs and Exp are each guaranteed to be between 0 and 1 (inclusive).
	// Therefore, the maximum Omega index is 1.0.
	// I think the minimum is -1.0, although not certain that can be attained.

	// First step is to identify all the distinct nodes
	set<NodeId> nodes;
	for(typename NodeToGroup<NodeId>::const_iterator i = ng1.begin(); i!=ng1.end(); i++) { nodes.insert(i->first); }
	for(typename NodeToGroup<NodeId>::const_iterator i = ng2.begin(); i!=ng2.end(); i++) { nodes.insert(i->first); }

	vector<NodeId> nodesv;
	for(typename set<NodeId>::const_iterator i=nodes.begin(); i!=nodes.end(); i++) { nodesv.push_back(*i); }
	const int N=nodesv.size();

	// Now to start populating the relevant statistics
	map<int,lli> A; // number of pairs of nodes sharing this number of clusters, where both covers agree on this number
	map<int,lli> N_bottom; // necessary for computing the expected number. For one cover, the distribution of the shared-cluster-count.
	map<int,lli> N_side;   // necessary for computing the expected number. For the other cover, the distribution of the shared-cluster-count.

	int minJK = 0;
	lli sumOfSquares = 0;
	for(int n=0; n<N; n++) {
		for(int m=0; m<n; m++) {
			// PP2(n,m);
			const NodeId n_ = nodesv.at(n);
			const NodeId m_ = nodesv.at(m);
			// PP2(n_,m_);
			const int a = ng1.sharedGroups(n_,m_);
			const int c = ng2.sharedGroups(n_,m_);
			// PP4(a,b,c,d);
			// assert(a!=1 || d!=0);
			if(a==c)
				A[a]++;
			N_bottom[a]++;
			N_side  [c]++;
			if(minJK < a)
				minJK=a;
			if(minJK < c)
				minJK=c;


			{ // Latouches's L2 norm
				// PP2_v(a-c, (a-c)*(a-c));
				sumOfSquares += (a-c)*(a-c);
			}
		}
	}
	PP1_v(minJK);

	lli bigN = 0;
	for(const pair<int,int> &Nj: N_bottom) {
		bigN += Nj.second;
	}
	{ // verification
		/*
		forEach(const typeof(pair<int,int>) &Aj, amd::mk_range(A)) {
			PP2(Aj.first, Aj.second);
		}
		forEach(const typeof(pair< pair<int,int>,int>) &Bij, amd::mk_range(B)) {
			PP3(Bij.first.first, Bij.first.second, Bij.second);
		}
		forEach(const typeof(pair<int,int>) &Nj, amd::mk_range(N_bottom)) {
			// PP2(Nj.first, Nj.second);
		}
		*/
		int verifyNumPairs3 = 0;
		for(const pair<int,int> &Nj: N_side  ) {
			// PP2(Nj.first, Nj.second);
			verifyNumPairs3 += Nj.second;
		}
		assert(bigN == N*(N-1)/2);
		assert(verifyNumPairs3 == N*(N-1)/2);
	}

	double L2norm = sqrt(sumOfSquares);

	lli numerator = 0;
	for(int j=0; j<=minJK; j++) {
		numerator += bigN * A[j];
		numerator -= N_bottom[j] * N_side[j];
	}
	lli denominator = 0;
	denominator += bigN * bigN;
	for(int j=0; j<=minJK; j++) {
		denominator -= N_bottom[j] * N_side[j];
		assert(denominator >= 0);
	}
	PP1_v(numerator);
	PP1_v(denominator);
	double O = double(numerator) / double(denominator);
	return make_pair(O,L2norm);
}

//! \brief Synchronize nodes of the destination collection with the nodes of
//! the source collection removing empty clusters and non-matching nodes
//!
//! \param dcls Grouping<NodeId>&  - reducing clusters
//! \param dnds unordered_set<NodeId>&  - reducing nodes
//! \param snds const unordered_set<NodeId>&  - node base
//! \return void
template <typename NodeId>
void syncNodes(Grouping<NodeId>& dcls, unordered_set<NodeId>& dnds, const unordered_set<NodeId>& snds)  // Sync nodes in the group 2
{
	dnds.clear();
	const auto sne = snds.end();
	// Start in the reverse direction for the vector to reduce the number of relocations
	for(auto icl = --dcls.end(); icl != --dcls.begin();) {
		for(auto ind = icl->begin(); ind != icl->end();) {
			// Remove all non-matching nodes
			if(snds.find(*ind) != sne)
				dnds.insert(*ind++);
			else ind = icl->erase(ind);
		}
		if(icl->empty())
			icl = --dcls.erase(icl);
		else --icl;
	}
}

template <typename NodeId>
void onmi(const char * file1, const char * file2, const char * syncfile
, const bool do_omega_also, const bool allnmis, float membership) {
	Grouping<NodeId>  g1, g2;
	if(syncfile) {
		unordered_set<NodeId>  nodes1;
		unordered_set<NodeId>  nodes2;
		g1 = fileToSet<NodeId>(file1, membership, &nodes1);
		g2 = fileToSet<NodeId>(file2, membership, &nodes2);
		if(nodes1.size() != nodes2.size()) {
            cerr << "WARNING, the number of nodes is different in the clusterings: "
                << nodes1.size() << " != " << nodes2.size() << ", synchronizing...\n";
			// Force sync with the first file wit the syncfile is not '-'
			if(nodes1.size() <= nodes2.size() || !(strlen(syncfile) == 1 && *syncfile == '-')) {
				syncNodes(g2, nodes2, nodes1);  // Sync nodes in the group 2 to nodes1 base
				//assert(nodes2.size() <= nodes1.size());
			} else {
				syncNodes(g1, nodes1, nodes2);  // Sync nodes in the group 1 to nodes2 base
				//assert(nodes1.size() <= nodes2.size());
			}
			//// Check whether the sync is successful
			//if(nodes1.size() != nodes2.size()) {
			//	cerr << "After the synchronization  ndsnum1: " << nodes1.size() << ", ndsnum2: " << nodes2.size() << endl;
			//	throw std::domain_error("Input clusterings have different node base and can't be synchronized gracefully\n");
			//}
		}
	} else {
		g1 = fileToSet<NodeId>(file1, membership);
		g2 = fileToSet<NodeId>(file2, membership);
	}

	PP1_v(g1.size());
	PP1_v(g2.size());
	unless(g1.size() > 0) throw EmptyFile();
	unless(g2.size() > 0) throw EmptyFile();
	NodeToGroup<NodeId> n2g1 = nodeToGroup(g1);
	NodeToGroup<NodeId> n2g2 = nodeToGroup(g2);
	if(n2g1.size() != n2g2.size())
		cerr << "WARNING, the number of nodes is different in the collections (the quality will be penalized): "
			<< n2g1.size() << " != " << n2g2.size() << endl;
	PP1_v(n2g1.size());
	PP1_v(n2g2.size());
	const OverlapMatrix om = overlapMatrix(n2g1, n2g2);

	OverlapMatrix omFlipped;
	omFlipped.N = om.N;
	for(pair< pair<int,int> ,int> p: om.om) {
		swap(p.first.first, p.first.second);
		bool wasInserted = omFlipped.om.insert(p).second;
		assert(wasInserted);
	}
	assert(omFlipped.om.size() == om.om.size());

	if(global_verbose_flag) {
		cout << "  \'" << file2 << "\' given \'" << file1 << "\"" << endl;
		for(int toId = 0; toId < (int)g2.size(); toId++) {
			PP(HX_given_BestY(omFlipped, g1, g2, toId));
		}
		cout << "  \'" << file1 << "\' given \'" << file2 << "\"" << endl;
		for(int fromId = 0; fromId < (int)g1.size(); fromId++) {
			PP(HX_given_BestY(om       , g2, g1, fromId));
		}
		cout << "Here:" << endl;
	}
	if(do_omega_also) {
		pair<double,double> OmegaAndL2Norm = omega(n2g1, n2g2);
		const double Omega = OmegaAndL2Norm.first;
		const double L2norm = OmegaAndL2Norm.second;
		cout << "Omega: " << Omega << " (L2norm: " << L2norm << "), ";
	}
	const auto  nmix = aaronNMI<Max>(om, omFlipped, g1, g2);  // NMImax
	if(nmix < numeric_limits<decltype(nmix)>::epsilon())
		throw domain_error("NMI is not applicable to the specified collections: 0, which says nothing about the similarity\n");
	if(allnmis) {
		cout << "NMImax: " << nmix
			<< ", NMIsqrt: " << aaronNMI<Sqrt>(om, omFlipped, g1, g2)
			<< ", NMIsum: " << aaronNMI<Sum>(om, omFlipped, g1, g2)
			<< ", NMIlfk: " << LFKNMI(om, omFlipped, g1, g2) << endl;
	} else {
		if(do_omega_also)
			cout << "NMImax: ";
		cout << nmix << endl;
	}
}
