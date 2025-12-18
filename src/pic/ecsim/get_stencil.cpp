
#include "indexer_4d.h"

#include <set>
#include <map>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdio>

#define _ECSIM_STENCIL_SELFTEST_

// ---- Self-check after indexer refactor --------------------------------------
// Compares the row assembled by GetStencil_import_cStencil(...) against the
// canonical ECSIM::GetStencil(...) implementation. Enabled when the macro
// _ECSIM_STENCIL_SELFTEST_ is defined (e.g., add -D_ECSIM_STENCIL_SELFTEST_).
// -----------------------------------------------------------------------------
#if defined(_ECSIM_STENCIL_SELFTEST_)
namespace {

  using LSCN =
    cLinearSystemCornerNode<PIC::Mesh::cDataCornerNode,
                            3, _PIC_STENCIL_NUMBER_, _PIC_STENCIL_NUMBER_+1, 16, 1, 1>;

  using RowEl  = LSCN::cMatrixRowNonZeroElementTable;
  using RhsEnt = LSCN::cRhsSupportTable;

  // Robust double comparison: relative to max(|a|,|b|) with an absolute floor
  inline bool near(double a, double b, double rel, double abs_floor=1e-300) {
    double scale = std::max(std::max(std::abs(a), std::abs(b)), abs_floor);
    return std::abs(a-b) <= rel * scale;
  }

  struct Key { // primary identity by logical address
    int i, j, k, c;    // c == iVar
    bool operator<(const Key& o) const {
      if (i!=o.i) return i<o.i;
      if (j!=o.j) return j<o.j;
      if (k!=o.k) return k<o.k;
      return c<o.c;
    }
  };

  // A single tap's full payload, aggregated if duplicates present
  struct AggTap {
    // numeric
    double value_sum = 0.0;      // sum of MatrixElementValue over duplicates
    int count = 0;

    // node identity (all duplicates must agree; we record the first)
    const void* node_ptr = nullptr;

    // parameter table aggregation/checks
    // for our case length=1, but keep generic
    std::vector<int> param_lengths;                  // per duplicate
    std::vector<std::vector<double>> param_values;   // per duplicate

    // support table (we only allow length<=1 today, but compare generically)
    std::multiset<const void*> support_ptrs; // multiset across duplicates
    std::vector<int> support_lengths;        // per duplicate
  };

  // Build map<Key, AggTap> from a flat row
  static std::map<Key, AggTap> index_flat_row(const RowEl* T, int n,
                                              /*collect params & supports*/ bool deep=true)
  {
    std::map<Key, AggTap> M;
    for (int t=0; t<n; ++t) {
      Key k{T[t].i, T[t].j, T[t].k, T[t].iVar};
      auto& a = M[k];
      a.value_sum += T[t].MatrixElementValue;
      a.count += 1;

      if (!a.node_ptr) a.node_ptr = (const void*)T[t].Node;

      if (deep) {
        a.param_lengths.push_back(T[t].MatrixElementParameterTableLength);
        if (T[t].MatrixElementParameterTableLength > 0) {
          std::vector<double> pv(T[t].MatrixElementParameterTableLength);
          for (int p=0; p<T[t].MatrixElementParameterTableLength; ++p)
            pv[p] = T[t].MatrixElementParameterTable[p];
          a.param_values.emplace_back(std::move(pv));
        } else {
          a.param_values.emplace_back();
        }

        a.support_lengths.push_back(T[t].MatrixElementSupportTableLength);
        for (int s=0; s<T[t].MatrixElementSupportTableLength; ++s)
          a.support_ptrs.insert(T[t].MatrixElementSupportTable[s]);
      }
    }
    return M;
  }

  struct MatrixDiff {
    std::vector<Key> onlyA;
    std::vector<Key> onlyB;
    struct BadVal { Key k; double a,b; } largest_val_delta;
    bool has_val_delta=false;

    std::vector<Key> node_mismatch;
    std::vector<Key> param_len_mismatch;
    std::vector<Key> param_val_mismatch;
    std::vector<Key> support_len_mismatch;
    std::vector<Key> support_ptr_mismatch;

    bool ok() const {
      return onlyA.empty() && onlyB.empty() && !has_val_delta &&
             node_mismatch.empty() && param_len_mismatch.empty() &&
             param_val_mismatch.empty() && support_len_mismatch.empty() &&
             support_ptr_mismatch.empty();
    }
  };

  // Compare two matrix rows with deep checks
  static MatrixDiff CompareMatrixEnhanced(const RowEl* A, int lenA,
                                          const RowEl* B, int lenB,
                                          double relTol_val,
                                          double relTol_param = 0.0) // params usually exact
  {
    MatrixDiff diff;
    auto MA = index_flat_row(A, lenA, /*deep=*/true);
    auto MB = index_flat_row(B, lenB, /*deep=*/true);

    // Discover missing/extra keys
    {
      auto ia=MA.begin(), ib=MB.begin();
      while (ia!=MA.end() || ib!=MB.end()) {
        if (ib==MB.end() || (ia!=MA.end() && ia->first<ib->first)) {
          diff.onlyA.push_back(ia->first); ++ia; continue;
        }
        if (ia==MA.end() || (ib!=MB.end() && ib->first<ia->first)) {
          diff.onlyB.push_back(ib->first); ++ib; continue;
        }
        // keys equal
        ++ia; ++ib;
      }
    }
    if (!diff.onlyA.empty() || !diff.onlyB.empty()) return diff;

    // Check payload for common keys
    for (auto& kv : MA) {
      const Key& k = kv.first;
      const AggTap& a = kv.second;
      const AggTap& b = MB[k];

      // value sum
      if (!near(a.value_sum, b.value_sum, relTol_val)) {
        double da = std::abs(a.value_sum - b.value_sum);
        if (!diff.has_val_delta || da > std::abs(diff.largest_val_delta.a - diff.largest_val_delta.b)) {
          diff.largest_val_delta = {k, a.value_sum, b.value_sum};
          diff.has_val_delta = true;
        }
      }

      // node pointer identity (all duplicates should correspond to same logical node)
      if (a.node_ptr != b.node_ptr) diff.node_mismatch.push_back(k);

      // parameter table length agreement per-duplicate (multiset compare)
      {
        auto sortv = [](std::vector<int>& v){ std::sort(v.begin(), v.end()); };
        auto LA = a.param_lengths, LB = b.param_lengths;
        sortv(LA); sortv(LB);
        if (LA != LB) diff.param_len_mismatch.push_back(k);
      }

      // parameter values equality (bag of vectors by length)
      {
        // multiset by (len, values)
        auto canon = [](const std::vector<std::vector<double>>& vv){
          std::vector<std::pair<int,std::vector<double>>> out;
          out.reserve(vv.size());
          for (auto& v : vv) out.emplace_back((int)v.size(), v);
          std::sort(out.begin(), out.end(),
                    [](auto& x, auto& y){
                      if (x.first!=y.first) return x.first<y.first;
                      return x.second<y.second;
                    });
          return out;
        };
        auto PA = canon(a.param_values);
        auto PB = canon(b.param_values);

        if (PA.size()!=PB.size()) diff.param_val_mismatch.push_back(k);
        else {
          for (size_t t=0; t<PA.size(); ++t) {
            if (PA[t].first != PB[t].first) { diff.param_val_mismatch.push_back(k); break; }
            int L = PA[t].first;
            for (int p=0; p<L; ++p) {
              if (!near(PA[t].second[p], PB[t].second[p], relTol_param)) {
                diff.param_val_mismatch.push_back(k); goto param_done;
              }
            }
          }
        }
      }
      param_done: ;

      // support length multiset equality
      {
        auto LA = a.support_lengths, LB = b.support_lengths;
        std::sort(LA.begin(), LA.end());
        std::sort(LB.begin(), LB.end());
        if (LA != LB) diff.support_len_mismatch.push_back(k);
      }

      // support pointer multiset equality
      if (a.support_ptrs != b.support_ptrs) diff.support_ptr_mismatch.push_back(k);
    }

    return diff;
  }

  struct RhsDiff {
    bool ptr_multiset_equal = true;
    bool coef_equal = true;
    // diagnostics
    std::map<const void*, int> onlyA_ptrs, onlyB_ptrs; // with multiplicities
    struct CoefDelta { const void* ptr; double a,b; };
    std::vector<CoefDelta> coef_deltas;
    bool ok() const { return ptr_multiset_equal && coef_equal; }
  };

  // Compare RHS (corner or center): pointers multiset + coefficients per pointer (sum)
  static RhsDiff CompareRhsEnhanced(const RhsEnt* A, int lenA,
                                    const RhsEnt* B, int lenB,
                                    double relTol_coef)
  {
    RhsDiff diff;

    // consolidate by pointer
    auto bucket = [](const RhsEnt* X, int n){
      std::map<const void*, std::pair<int,double>> M; // ptr -> (count, sumCoef)
      for (int t=0; t<n; ++t) {
        auto& e = M[X[t].AssociatedDataPointer];
        e.first += 1;
        e.second += X[t].Coefficient;
      }
      return M;
    };

    auto MA = bucket(A, lenA);
    auto MB = bucket(B, lenB);

    // pointer multiset equality (counts)
    // also track missing/extra
    for (auto& kv : MA) {
      auto it = MB.find(kv.first);
      if (it==MB.end()) { diff.ptr_multiset_equal = false; diff.onlyA_ptrs[kv.first] = kv.second.first; }
      else if (it->second.first != kv.second.first) { diff.ptr_multiset_equal = false; }
    }
    for (auto& kv : MB) {
      auto it = MA.find(kv.first);
      if (it==MA.end()) { diff.ptr_multiset_equal = false; diff.onlyB_ptrs[kv.first] = kv.second.first; }
    }

    // coefficients per pointer (sum)
    for (auto& kv : MA) {
      auto it = MB.find(kv.first);
      if (it==MB.end()) continue;
      double a = kv.second.second, b = it->second.second;
      if (!near(a,b,relTol_coef)) {
        diff.coef_equal = false;
        diff.coef_deltas.push_back({kv.first,a,b});
      }
    }

    return diff;
  }

  // Pretty-printers for diagnostics
  static void print_key(const char* prefix, const Key& k) {
    fprintf(stderr, "%s(i,j,k,c)=(%d,%d,%d,%d)\n", prefix, k.i,k.j,k.k,k.c);
  }

  static void report_matrix_diff(const MatrixDiff& D) {
    if (!D.onlyA.empty() || !D.onlyB.empty()) {
      if (!D.onlyA.empty()) {
        fprintf(stderr,"    Matrix: keys only in ROW (missing in REF): %zu (show up to 12)\n", D.onlyA.size());
        for (size_t t=0;t<std::min<size_t>(12,D.onlyA.size());++t) print_key("      ", D.onlyA[t]);
      }
      if (!D.onlyB.empty()) {
        fprintf(stderr,"    Matrix: keys only in REF (missing in ROW): %zu (show up to 12)\n", D.onlyB.size());
        for (size_t t=0;t<std::min<size_t>(12,D.onlyB.size());++t) print_key("      ", D.onlyB[t]);
      }
    }
    if (D.has_val_delta) {
      auto& dv = D.largest_val_delta;
      fprintf(stderr,"    Matrix: largest value delta at "); print_key("", dv.k);
      fprintf(stderr,"      row=% .6e  ref=% .6e  |Δ|=% .6e\n",
              dv.a, dv.b, std::abs(dv.a-dv.b));
    }
    auto dump_keys = [](const char* hdr, const std::vector<Key>& V){
      if (V.empty()) return;
      fprintf(stderr,"    %s: %zu (show up to 12)\n", hdr, V.size());
      for (size_t t=0;t<std::min<size_t>(12,V.size());++t) print_key("      ", V[t]);
    };
    dump_keys("Matrix: node pointer mismatch", D.node_mismatch);
    dump_keys("Matrix: parameter LENGTH mismatch", D.param_len_mismatch);
    dump_keys("Matrix: parameter VALUE mismatch", D.param_val_mismatch);
    dump_keys("Matrix: support LENGTH mismatch", D.support_len_mismatch);
    dump_keys("Matrix: support POINTER mismatch", D.support_ptr_mismatch);
  }

  static void report_rhs_diff(const char* tag, const RhsDiff& D) {
    if (!D.ptr_multiset_equal) {
      fprintf(stderr, "    %s: pointer multiset mismatch\n", tag);
      if (!D.onlyA_ptrs.empty()) {
        fprintf(stderr, "      only in ROW (ptr,count): show up to 12\n");
        int n=0; for (auto& kv : D.onlyA_ptrs) { if (n++>=12) break; fprintf(stderr,"        %p x%d\n", kv.first, kv.second); }
      }
      if (!D.onlyB_ptrs.empty()) {
        fprintf(stderr, "      only in REF (ptr,count): show up to 12\n");
        int n=0; for (auto& kv : D.onlyB_ptrs) { if (n++>=12) break; fprintf(stderr,"        %p x%d\n", kv.first, kv.second); }
      }
    }
    if (!D.coef_equal) {
      fprintf(stderr, "    %s: coefficient mismatch for %zu pointers (show up to 12)\n",
              tag, D.coef_deltas.size());
      for (size_t t=0;t<std::min<size_t>(12, D.coef_deltas.size()); ++t) {
        auto& cd = D.coef_deltas[t];
        fprintf(stderr, "        ptr=%p  row=% .6e  ref=% .6e  |Δ|=% .6e\n",
                cd.ptr, cd.a, cd.b, std::abs(cd.a-cd.b));
      }
    }
  }

  // ----- SELF-TEST DRIVER -------------------------------------------------------
static void PostAssembleSelfCheck(
    int i, int j, int k, int iVar,
    RowEl* rowSet, int rowNze, double rowRhs,
    RhsEnt* rowRhsCorner, int rowRhsCornerLen,
    RhsEnt* rowRhsCenter, int rowRhsCenterLen,
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node)
{
  // Build reference via canonical assembler
  RowEl  refSet[_PIC_STENCIL_NUMBER_];      int refNze = 0;     double refRhs = 0.0;
  RhsEnt refCorner[_PIC_STENCIL_NUMBER_+1]; int refCornerLen = 0;
  RhsEnt refCenter[16];                     int refCenterLen = 0;

  vector<RhsEnt> support_corner_vector,support_center_vector; 

  PIC::FieldSolver::Electromagnetic::ECSIM::GetStencil(
      i, j, k, iVar,
      refSet, refNze, refRhs,
      refCorner, refCornerLen,
      refCenter, refCenterLen,
      node,support_corner_vector,support_center_vector);

  CompareArray_Unordered<RhsEnt>(refCorner,refCornerLen,rowRhsCorner,rowRhsCornerLen,true,"refCorner"); 
  CompareArray_Unordered<RhsEnt>(refCenter,refCenterLen,rowRhsCenter,rowRhsCenterLen,true,"refCorner");

  //RhsEnt::CompareArray(refCorner,rowRhsCorner,rowRhsCornerLen,true); 
  
  //RowEl::CompareArray(refSet,rowSet,rowNze,true);
   
  CompareArray_Unordered<RowEl>(rowSet,rowNze,refSet,refNze,true,"sdf");

  //CompareArray_Unordered<RhsEnt>(rowRhsCorner,rowRhsCornerLen,refCorner,refCornerLen,true,"RHS unmatched");
  //CompareArray_Unordered<RhsEnt>(rowRhsCenter,rowRhsCenterLen,refCenter,refCenterLen,true,"RHS unmatched"); 

  // Tolerances
  const double relTol_val   = 1e-12;  // matrix numeric
  const double relTol_param = 0.0;    // parameter tables (expect exact)
  const double relTol_rhs   = 1e-14;  // RHS coeff sums per pointer
  const double relTol_rhs_s = 1e-14;  // scalar RHS

  // Deep comparisons (matrix + RHS)
  auto mdiff = CompareMatrixEnhanced(rowSet, rowNze,
                                     refSet,  refNze,
                                     relTol_val, relTol_param);

  auto cdiff = CompareRhsEnhanced(rowRhsCorner, rowRhsCornerLen,
                                  refCorner,    refCornerLen,
                                  relTol_rhs);

  auto zdiff = CompareRhsEnhanced(rowRhsCenter, rowRhsCenterLen,
                                  refCenter,     refCenterLen,
                                  relTol_rhs);

  bool okRhsScalar = near(rowRhs, refRhs, relTol_rhs_s, 1e-300);

  if (!(mdiff.ok() && cdiff.ok() && zdiff.ok() && okRhsScalar)) {
    std::fprintf(stderr,
      "\n[ECSIM::StencilSelfTest] Mismatch at (i,j,k,c)=(%d,%d,%d,%d)\n"
      "  row:    nze=%d  rhs=%.16e   rhsCorner=%d  rhsCenter=%d\n"
      "  refrow: nze=%d  rhs=%.16e   rhsCorner=%d  rhsCenter=%d\n",
      i, j, k, iVar,
      rowNze, rowRhs, rowRhsCornerLen, rowRhsCenterLen,
      refNze, refRhs, refCornerLen, refCenterLen);

    if (!mdiff.ok())  report_matrix_diff(mdiff);
    if (!cdiff.ok())  report_rhs_diff("RHS Corner", cdiff);
    if (!zdiff.ok())  report_rhs_diff("RHS Center", zdiff);
    if (!okRhsScalar) {
      std::fprintf(stderr, "    RHS scalar mismatch: row=%.16e  ref=%.16e  |Δ|=%.3e\n",
                   rowRhs, refRhs, std::abs(rowRhs - refRhs));
    }

    // Project's fatal exit (same pattern used elsewhere in AMPS)
    exit(__LINE__, __FILE__, "ECSIM stencil self-test failed");
  }
}


} // namespace
#endif // _ECSIM_STENCIL_SELFTEST_


void PIC::FieldSolver::Electromagnetic::ECSIM::GetStencil_import_cStencil(int i,int j,int k,int iVar,cLinearSystemCornerNode<PIC::Mesh::cDataCornerNode,3,_PIC_STENCIL_NUMBER_,_PIC_STENCIL_NUMBER_+1,16,1,1>::cMatrixRowNonZeroElementTable* MatrixRowNonZeroElementTable,int& NonZeroElementsFound,double& rhs,
			     cLinearSystemCornerNode<PIC::Mesh::cDataCornerNode,3,_PIC_STENCIL_NUMBER_,_PIC_STENCIL_NUMBER_+1,16,1,1>::cRhsSupportTable* RhsSupportTable_CornerNodes,int& RhsSupportLength_CornerNodes,
			     cLinearSystemCornerNode<PIC::Mesh::cDataCornerNode,3,_PIC_STENCIL_NUMBER_,_PIC_STENCIL_NUMBER_+1,16,1,1>::cRhsSupportTable* RhsSupportTable_CenterNodes,int& RhsSupportLength_CenterNodes, 
			     cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node,                                                                           vector<cLinearSystemCornerNode<PIC::Mesh::cDataCornerNode,3,_PIC_STENCIL_NUMBER_,_PIC_STENCIL_NUMBER_+1,16,1,1>::cRhsSupportTable>& support_corner_vector,
              vector<cLinearSystemCornerNode<PIC::Mesh::cDataCornerNode,3,_PIC_STENCIL_NUMBER_,_PIC_STENCIL_NUMBER_+1,16,1,1>::cRhsSupportTable>& support_center_vector) {
  
  using namespace PIC::FieldSolver::Electromagnetic::ECSIM;    


//  GetStencil(i,j,k,iVar,MatrixRowNonZeroElementTable,NonZeroElementsFound,rhs,RhsSupportTable_CornerNodes,RhsSupportLength_CornerNodes,RhsSupportTable_CenterNodes,RhsSupportLength_CenterNodes,node);
//  return;
  
  // No.0-No.26  stencil Ex
  // No.27-No.53 stencil Ey
  // No.54-No.80 stencil Ez
  
  // i+indexadd[ii](ii:0,1,2), j+indexAdd[jj](jj:0,1,2), k+indexAdd[kk](kk:0,1,2)
  // index number: ii+3*jj+9*kk 
  // No.0: i,j,k No.1 i-1,j,k No.2 i+1,j,k
  // No.3: i,j-1,k No.4:i-1,j-1,k No.5 i+1,j-1,k 
 
  // double cLighgt, dt;
  double x[3];
  //power of 3 array created
  //for some pgi compilers that cannot convert result of pow() from double to int correctly
  int pow3_arr[3]={1,3,9};
  int index[3] = {i,j,k};
  int nCell[3] = {_BLOCK_CELLS_X_,_BLOCK_CELLS_Y_,_BLOCK_CELLS_Z_};
  double dx[3],coeff[3],coeffSqr[3],coeff4[3]; 

  //init indexers 
  //static 
	  cIndexer4D indexer[3]; //iVar=0...2 

//  static cIndexer4D indexer_rhs;

  if (indexer[0].isInitialized()==false) {
    //init indexer

    int min_limit[4]={-5,-5,-5,0},max_limit[4]={5,5,5,2};

    for (int i=0;i<3;i++) indexer[i].init(min_limit,max_limit);

 //   indexer_rhs.init(min_limit,max_limit);
  }

for (int i=0;i<3;i++) indexer[i].reset();

  // Canonical mapping from physical offsets (di,dj,dk) and unknown component
  // (iVarIndex: 0->Ex, 1->Ey, 2->Ez) into the compact element slot *for this row iVar*.
  // This keeps iElement <-> (di,dj,dk,iVarIndex) consistent across all terms.
  auto map_iElement = [&](int di, int dj, int dk, int iVarIndex)->int {
    return indexer[iVar].get_idx(di, dj, dk, iVarIndex);
  };

  //------------------------------------------------------------------------------
  // ensure_curlcurl_stencil_initialized
  //------------------------------------------------------------------------------
  // PURPOSE
  //   One-shot initialization wrapper for the 4th-order curl–curl(E) stencil.
  //   On the *first* call, it invokes the provided initializer
  //
  //     PIC::FieldSolver::Electromagnetic::ECSIM::Stencil::FourthOrder::
  //     InitCurlCurlEStencils(cCurlCurlEStencil* CurlCurl, double dx, double dy, double dz)
  //
  //   with the given spacings (dx,dy,dz) and caches the result in a static
  //   object. On subsequent calls, it *does not* re-initialize; it simply copies
  //   the cached stencil into the caller-provided `out`.
  //
  // NOTES
  //   - This is intentionally minimal: a single static boolean gate. If your
  //     grid spacings can change at runtime, extend this to also cache (dx,dy,dz)
  //     and rebuild when they differ. As written, it initializes once and reuses.
  //   - Copying `s_stencil` into `out` relies on `cStencil`/`cCurlCurlEStencil`
  //     being trivially copyable or having proper copy semantics.
  //   - Not thread-safe. If used from multiple threads, guard the first-call
  //     section with a mutex.
  //
  // Example usage:
  // cCurlCurlEStencil *curlcurl;
  // curlcurl=get_curlcurl_stencil(dx, dy, dz);
  //------------------------------------------------------------------------------
  auto get_grad_div_stencil =
    [&]() -> PIC::FieldSolver::Electromagnetic::ECSIM::Stencil::cGradDivEStencil* {
      using PIC::FieldSolver::Electromagnetic::ECSIM::Stencil::FourthOrder::InitGradDivEBStencils;

      static bool               s_inited   = false;
      static PIC::FieldSolver::Electromagnetic::ECSIM::Stencil::cGradDivEStencil s_stencil[3];          // cached once

      if (!s_inited) {
        InitGradDivEBStencils(s_stencil, 1.0,1.0,1.0);  // build exactly once
        s_inited = true;
      }

      // hand back the cached stencil to the caller
      return s_stencil;
    };

  auto get_laplacian_stencil =
    [&]() -> PIC::FieldSolver::Electromagnetic::ECSIM::Stencil::cLaplacianStencil* {
      using PIC::FieldSolver::Electromagnetic::ECSIM::Stencil::FourthOrder::InitLaplacianStencil;

      static bool               s_inited   = false;
      static PIC::FieldSolver::Electromagnetic::ECSIM::Stencil::cLaplacianStencil  s_stencil;          // cached once

      if (!s_inited) {
        InitLaplacianStencil(&s_stencil, 1.0,1.0,1.0);  // build exactly once
        s_inited = true;
      }

      // hand back the cached stencil to the caller
      return &s_stencil;
    };

  //------------------------------------------------------------------------------
  // init_metrics_and_time_factors
  //------------------------------------------------------------------------------
  // PURPOSE
  //   Compute per-dimension grid metrics and ECSIM time-centering factors used
  //   throughout stencil assembly:
  //     • dx[d]      : physical cell size along axis d (after unit conversion)
  //     • x[d]       : coordinate of the current corner index along axis d
  //     • coeff[d]   : θ c Δt / dx[d]           (single-derivative scale factor)
  //     • coeffSqr[d]: (θ c Δt / dx[d])^2       (second-derivative scale factor)
  //     • coeff4[d]  : (θ c Δt / dx[d]) / 4     (used for face-avg in curl(B))
  //
  // INPUTS (must exist in the enclosing scope):
  //   node   → current AMR tree node (provides xmin/xmax bounds)
  //   nCell  → {NX, NY, NZ} number of cells per block in each dimension
  //   index  → {i, j, k} integer corner indices within the block
  //   cDt    → c * Δt  (speed of light times timestep)
  //   theta  → ECSIM time-centering parameter
  //   length_conv → optional unit conversion factor for lengths
  //
  // OUTPUTS (filled by this lambda; exist in enclosing scope):
  //   dx[3], x[3], coeff[3], coeffSqr[3], coeff4[3]
  //
  // NOTES
  //   • x[d] is evaluated at the *corner* location corresponding to (i,j,k).
  //   • coeff and coeffSqr set the correct scaling for first/second derivatives
  //     in discrete operators (e.g., grad, div, Laplacian, grad–div).
  //   • coeff4 is exactly the 1/4 factor used when averaging 4 center samples
  //     on a face to approximate curl(B) with a compact stencil.
  //------------------------------------------------------------------------------
  auto init_metrics_and_time_factors = [&]() {
    for (int iDim = 0; iDim < 3; ++iDim) {
      // Raw cell size in node coordinates
      dx[iDim] = (node->xmax[iDim] - node->xmin[iDim]) / nCell[iDim];

      // Apply unit conversion (if any)
      dx[iDim] *= length_conv;

      // Corner coordinate along this axis at index[iDim]
      x[iDim] = node->xmin[iDim]
            + index[iDim] * (node->xmax[iDim] - node->xmin[iDim]) / nCell[iDim];

      // Time-centering / metric factors:
      //   coeff    = θ c Δt / Δx
      //   coeffSqr = (θ c Δt / Δx)^2
      //   coeff4   = (θ c Δt / Δx) / 4 (face-average factor for curl(B))
      coeff[iDim]    = (cDt / dx[iDim]) * theta;
      coeffSqr[iDim] = coeff[iDim] * coeff[iDim];
      coeff4[iDim]   = 0.25 * coeff[iDim];
    }
  };

  init_metrics_and_time_factors();


  //------------------------------------------------------------------------------
  // boundary_short_circuit
  //------------------------------------------------------------------------------
  // PURPOSE
  //   Handle non-periodic physical boundaries for the ΔE solve at corner (i,j,k).
  //   If (i,j,k) is a boundary corner, we enforce the Dirichlet condition ΔE=0 by
  //   emitting a *single diagonal row* and skipping all other assembly. This
  //   makes the row:
  //        A_00 := 1  (stored via PARAMETER[0] = 1.0; numeric value filled later)
  //        rhs  := 0
  //   and sets lengths of RHS support tables to zero.
  //
  // WHAT IT WRITES (index 0 only):
  //   MatrixRowNonZeroElementTable[0].{i,j,k}            ← (i,j,k)
  //   MatrixRowNonZeroElementTable[0].iVar               ← iVar (row/equation component)
  //   MatrixRowNonZeroElementTable[0].MatrixElementValue ← 0.0      (value is filled later)
  //   MatrixRowNonZeroElementTable[0].MatrixElementParameterTable[0] ← 1.0 (unit diagonal)
  //   MatrixRowNonZeroElementTable[0].MatrixElementParameterTableLength ← 1
  //   MatrixRowNonZeroElementTable[0].MatrixElementSupportTableLength   ← 0  (no supports)
  //   MatrixRowNonZeroElementTable[0].MatrixElementSupportTable[0]      ← NULL
  //   MatrixRowNonZeroElementTable[0].Node ← node
  //
  // SIDE EFFECTS:
  //   NonZeroElementsFound      ← 1 (or 0 for the "right" boundary corner case)
  //   rhs                       ← 0
  //   RhsSupportLength_CenterNodes ← 0
  //   RhsSupportLength_CornerNodes ← 0
  //
  // RETURN VALUE:
  //   true  → boundary handled (caller should `return;` from GetStencil)
  //   false → not a boundary (caller continues with normal assembly)
  //
  // NOTES:
  //   • Active only when periodic BCs are OFF: _PIC_BC__PERIODIC_MODE_==_PIC_BC__PERIODIC_MODE_OFF_
  //   • The special case `isRightBoundaryCorner(x,node)` mirrors the original code’s
  //     treatment where the single entry is suppressed by setting NonZeroElementsFound=0.
  //------------------------------------------------------------------------------
  auto boundary_short_circuit = [&]() -> bool {
    if (_PIC_BC__PERIODIC_MODE_ == _PIC_BC__PERIODIC_MODE_OFF_)
      if (isBoundaryCorner(x, node)) {
        MatrixRowNonZeroElementTable[0].i = i;
        MatrixRowNonZeroElementTable[0].j = j;
        MatrixRowNonZeroElementTable[0].k = k;

        MatrixRowNonZeroElementTable[0].MatrixElementValue = 0.0; // numeric set later by updater
        MatrixRowNonZeroElementTable[0].iVar = iVar;

        // One-parameter design: put 1.0 on PARAMETER[0] to encode the unit diagonal.
        MatrixRowNonZeroElementTable[0].MatrixElementParameterTable[0] = 1.0;
        MatrixRowNonZeroElementTable[0].MatrixElementParameterTableLength = 1;

        // No supports for boundary diagonal.
        MatrixRowNonZeroElementTable[0].MatrixElementSupportTableLength = 0;
        MatrixRowNonZeroElementTable[0].MatrixElementSupportTable[0] = NULL;

        MatrixRowNonZeroElementTable[0].Node = node;

        // Set row sizes / rhs
        NonZeroElementsFound = 1;
        rhs = 0.0;
        RhsSupportLength_CenterNodes = 0;
        RhsSupportLength_CornerNodes = 0;

        // Preserve original behavior: suppress even this entry on the "right" boundary.
        if (isRightBoundaryCorner(x, node)) NonZeroElementsFound = 0;

        return true; // caller should `return;`
      }
    return false; // not handled; proceed with normal stencil assembly
  };

 if (boundary_short_circuit()) return;



  if (!initMassMatrixOffsetTable) computeMassMatrixOffsetTable(); 

  const int indexAddition[3] = {0,-1,1};
  const int reversed_indexAddition[3] = {1,0,2};  //table to determine ii,jj,kk from i,j,k 

  char * NodeDataOffset = node->block->GetCornerNode(_getCornerNodeLocalNumber(i,j,k))->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset;

  /*
  emit_compact_mass_structure
  -----------------------------------------
  PURPOSE
    Populate the 3×3×3 *compact* (radius-1) coupling structure of the ECSIM
    ΔE system’s **mass matrix term** for a single row (equation component p=iVar)
    at corner (i,j,k). For each target/column component q∈{Ex,Ey,Ez} and each
    neighbor corner (i+αx, j+αy, k+αz), α•∈{−1,0,+1}, the lambda:
      1) records WHICH column unknown it couples to (indices + component),
      2) allocates one PARAMETER slot (kept at length=1, value set later),
      3) stores one SUPPORT pointer to the mass/metric entry M_{pq} for this tap.

  MATHEMATICAL ROLE
    This is the structural encoding of the compact product (per row component p):
      \[
        (M\,\Delta\mathbf{E})_p(i,j,k)
        = \sum_{q\in\{x,y,z\}}
          \sum_{\alpha_x,\alpha_y,\alpha_z\in\{-1,0,1\}}
          M_{pq}(i,j,k;\boldsymbol{\alpha})\;
          \Delta E_q(i+\alpha_x,\;j+\alpha_y,\;k+\alpha_z).
      \]
    In the full operator
      \[
        \big[I + (\theta c\Delta t)^2(\nabla\nabla\!\cdot - \nabla^2 I)
          + 4\pi\,\theta\,\Delta t\,M\big]\;\Delta\mathbf{E}=\mathbf{r},
      \]
    this mass product is later scaled by \(4\pi\,\theta\,\Delta t\) and added to
    the numeric matrix entry during the **matrix-update** pass.

  WHAT GETS WRITTEN PER NONZERO (iElement)
    • Column addressing:
        .i,.j,.k   ← (i+αx, j+αy, k+αz),
        .iVar      ← q (the column component).
    • Value placeholder:
        .MatrixElementValue ← 0.0 (the actual number is computed later).
    • PARAMETER (for I + curl–curl part; original design = one scalar):
        .MatrixElementParameterTable[0] ← 0.0;
        .MatrixElementParameterTableLength ← 1;
      The matrix-update routine will overwrite/accumulate this to form the
      compact \(I + (\nabla\nabla\!\cdot - \nabla^2 I)\) contribution.
    • SUPPORT (for mass matrix M; original design = one pointer):
        .MatrixElementSupportTableLength ← 1;
        .MatrixElementSupportTable[0] ←
            (double*)NodeDataOffset + MassMatrixOffsetIndex
          + MassMatrixOffsetTable[iVar][iElement];
      This is an opaque address to the stored coefficient
      \(M_{pq}(i,j,k;\boldsymbol{\alpha})\). The update pass dereferences it and
      adds \( (4\pi\,\theta\,\Delta t)\times M_{pq} \) to the numeric entry.

  LAYOUT / ORDERING
    • 27 neighbors per component, 3 components → **81** entries per interior row.
      Linear index: iElement = q*27 + (ii + 3*jj + 9*kk), with ii,jj,kk∈{0,1,2}.
    • (p=iVar) identifies the row component; (q=iVarIndex) identifies the column
      component block. The per-tap mass offset uses both: MassMatrixOffsetTable[p][iElement].

  NOTES / GOTCHAS
    • Do not assign numeric coefficients here—leave MatrixElementValue at 0.0.
      Numbers are produced by the matrix-update routine using PARAMETER[0] & SUPPORT[0].
    • Pointer arithmetic assumes NodeDataOffset is aligned to “double*” slots for this
      region. If your build uses byte offsets, switch to `char*` arithmetic consistently.
    • Boundary/periodic remapping and row compaction are performed elsewhere.
*/
  auto emit_compact_mass_structure = [&]() {
    for (int iVarIndex = 0; iVarIndex < 3; iVarIndex++) {
      for (int ii = 0; ii < 3; ii++) {
        for (int jj = 0; jj < 3; jj++) {
          for (int kk = 0; kk < 3; kk++) {
            int iNode = i + indexAddition[ii];
            int jNode = j + indexAddition[jj];
            int kNode = k + indexAddition[kk];

            // Linear slot: 27 per component, 3 components total → 81 entries
            int iElement = iVarIndex * 27 + ii + jj * 3 + kk * 9;
	    int legacyElem = iVarIndex * 27 + ii + jj * 3 + kk * 9;

            const int di = indexAddition[ii];
            const int dj = indexAddition[jj];
            const int dk = indexAddition[kk];
            indexer[iVar].check_and_set(iElement,di, dj, dk, iVarIndex);
	    iElement=indexer[iVar].get_idx(di, dj, dk, iVarIndex); 


            MatrixRowNonZeroElementTable[iElement].i = iNode;
            MatrixRowNonZeroElementTable[iElement].j = jNode;
            MatrixRowNonZeroElementTable[iElement].k = kNode;

            // already initialized in LinearSystemSolver->h
            MatrixRowNonZeroElementTable[iElement].MatrixElementValue = 0.0;
            MatrixRowNonZeroElementTable[iElement].iVar = iVarIndex;

            // Original design: exactly one parameter per element
            MatrixRowNonZeroElementTable[iElement].MatrixElementParameterTable[0] = 0.0;
            MatrixRowNonZeroElementTable[iElement].MatrixElementParameterTableLength = 1;

            // Exactly one support pointer per element (mass/metric slot)
            MatrixRowNonZeroElementTable[iElement].MatrixElementSupportTableLength = 1;

            // Pointer math as in the original codebase:
            MatrixRowNonZeroElementTable[iElement].MatrixElementSupportTable[0] =
              (double*)NodeDataOffset
              + MassMatrixOffsetIndex
              + MassMatrixOffsetTable[iVar][legacyElem];
          }
        }
      }
    }
  };

  // invoke it
  emit_compact_mass_structure();
    

  int indexOffset[5] = {0,-1,1,-2,2};

  int reversed_indexOffset[5]={3,1,0,2,4}; //table to convert i,j,k ->ii,jj,kk

if (_PIC_STENCIL_NUMBER_==375) {
  for (int iVarIndex=0; iVarIndex<3; iVarIndex++){
    int cntTemp =0;
    
    for (int kk=0;kk<5;kk++){
      for (int jj=0;jj<5;jj++){	
	for (int ii=0;ii<5;ii++){       
	  
	  if (ii<3 && jj<3 && kk<3) continue;
          int iNode = i+indexOffset[ii];
          int jNode = j+indexOffset[jj];
          int kNode = k+indexOffset[kk];
          int iElement = 81+iVarIndex*98+cntTemp;

      {
        const int di = indexOffset[ii];
        const int dj = indexOffset[jj];
        const int dk = indexOffset[kk];
        //indexer[iVar].check_and_set(iElement,di, dj, dk, iVarIndex);
	iElement=indexer[iVar].get_idx(di, dj, dk, iVarIndex);
      }

	  cntTemp++;
          MatrixRowNonZeroElementTable[iElement].i=iNode;
          MatrixRowNonZeroElementTable[iElement].j=jNode;
          MatrixRowNonZeroElementTable[iElement].k=kNode;

          MatrixRowNonZeroElementTable[iElement].MatrixElementValue=0.0;
          MatrixRowNonZeroElementTable[iElement].iVar=iVarIndex;
          MatrixRowNonZeroElementTable[iElement].MatrixElementParameterTable[0]=0.0;
          MatrixRowNonZeroElementTable[iElement].MatrixElementParameterTableLength=1;
          MatrixRowNonZeroElementTable[iElement].MatrixElementSupportTableLength = 0;
        }
      }
    }
  }
}


  //laplacian
  auto *laplacian=get_laplacian_stencil();

  for (int idim=0;idim<3;idim++) {
    cStencil::cStencilData *st=LaplacianStencil+idim;

    for (int it=0;it<st->Length;it++) {
      int ii=reversed_indexAddition[st->Data[it].i+1]; 
      int jj=reversed_indexAddition[st->Data[it].j+1];
      int kk=reversed_indexAddition[st->Data[it].k+1];

      int index=ii+jj*3+kk*9;
      int iElement = index + iVar*27;

      {
        const int di = st->Data[it].i;
        const int dj = st->Data[it].j;
        const int dk = st->Data[it].k;
        //indexer[iVar].check_and_set(iElement,di, dj, dk, iVar); 
	iElement=indexer[iVar].get_idx(di, dj, dk, iVar);
      }

      //minus laplacian
      MatrixRowNonZeroElementTable[iElement].MatrixElementParameterTable[0]-=st->Data[it].a*coeffSqr[idim]; 
    }
  }
  

  //plus self
  //
  
int idx1;
      {
        const int di = 0;
        const int dj = 0;
        const int dk = 0;
        //indexer[iVar].check_and_set(27*iVar,di, dj, dk, iVar);
	idx1=indexer[iVar].get_idx(di, dj, dk, iVar);
      }
  
  MatrixRowNonZeroElementTable[idx1].MatrixElementParameterTable[0]+=1;


  //plus graddiv E
  for (int iVarIndex=0;iVarIndex<3;iVarIndex++) {
    cStencil::cStencilData *st=&GradDivStencil[iVar][iVarIndex];
    cStencil::cStencilData *st375=&GradDivStencil375[iVar][iVarIndex];

    for (int it=0;it<st->Length;it++) {
      int ii=reversed_indexAddition[st->Data[it].i+1];
      int jj=reversed_indexAddition[st->Data[it].j+1];
      int kk=reversed_indexAddition[st->Data[it].k+1];

      int nodeIndex=ii+jj*3+kk*9;
      int iElement = nodeIndex + iVarIndex*27;

      {
        const int di = st->Data[it].i;
        const int dj = st->Data[it].j;
        const int dk = st->Data[it].k;
        //indexer[iVar].check_and_set(iElement,di, dj, dk, iVarIndex);
	iElement=indexer[iVar].get_idx(di, dj, dk, iVarIndex);
      }

      MatrixRowNonZeroElementTable[iElement].MatrixElementParameterTable[0]+=(1-corrCoeff)*st->Data[it].a*coeff[iVar]*coeff[iVarIndex];
    }


    //add contribution from a larger stencil (divE correction)
    for (int it=0;it<st375->Length;it++) {
      if ((st375->Data[it].i<-1)||(st375->Data[it].i>1) || (st375->Data[it].j<-1)||(st375->Data[it].j>1) ||(st375->Data[it].k<-1)||(st375->Data[it].k>1) ) continue;

      int ii=reversed_indexAddition[st375->Data[it].i+1];
      int jj=reversed_indexAddition[st375->Data[it].j+1];
      int kk=reversed_indexAddition[st375->Data[it].k+1];

      int nodeIndex=ii+jj*3+kk*9;
      int iElement = nodeIndex + iVarIndex*27;

      {
        const int di = st375->Data[it].i;
        const int dj = st375->Data[it].j;
        const int dk = st375->Data[it].k;
        //indexer[iVar].check_and_set(iElement,di, dj, dk, iVarIndex);
	iElement=indexer[iVar].get_idx(di, dj, dk, iVarIndex);
      }

      MatrixRowNonZeroElementTable[iElement].MatrixElementParameterTable[0]+=corrCoeff*st375->Data[it].a*coeff[iVar]*coeff[iVarIndex];
    }
  }


int OrderingOffsetTable[5][5][5];

if (_PIC_STENCIL_NUMBER_==375) {   
  int cntTemp = 0;

  for (int kk=0;kk<5;kk++){
    for (int jj=0;jj<5;jj++){
      for (int ii=0;ii<5;ii++){
        if (ii<3 && jj<3 && kk<3) continue;

        OrderingOffsetTable[ii][jj][kk]=cntTemp;
        cntTemp++;
      }
    }
  }

  for (int iVarIndex=0;iVarIndex<3;iVarIndex++){
    cStencil::cStencilData *st375=&GradDivStencil375[iVar][iVarIndex];

    for (int it=0;it<st375->Length;it++) {

      int ii=reversed_indexOffset[st375->Data[it].i+2];
      int jj=reversed_indexOffset[st375->Data[it].j+2];
      int kk=reversed_indexOffset[st375->Data[it].k+2];

      if (ii<3 && jj<3 && kk<3) continue;
      int iElement = 81+iVarIndex*98+OrderingOffsetTable[ii][jj][kk];
      int nodeIndex= 27+OrderingOffsetTable[ii][jj][kk];

      {
        const int di = st375->Data[it].i;
        const int dj = st375->Data[it].j;
        const int dk = st375->Data[it].k;
        //indexer[iVar].check_and_set(iElement,di, dj, dk, iVarIndex);
	iElement=indexer[iVar].get_idx(di, dj, dk, iVarIndex);
      }

      MatrixRowNonZeroElementTable[iElement].MatrixElementParameterTable[0]+=corrCoeff*st375->Data[it].a*coeff[iVar]*coeff[iVarIndex];
    }
  }
}


  //find corners outside the boundary
  vector<int> pointLeft;
  int kMax=_BLOCK_CELLS_Z_,jMax=_BLOCK_CELLS_Y_,iMax=_BLOCK_CELLS_X_;

  for (int ii=0;ii<_PIC_STENCIL_NUMBER_;ii++) {
    MatrixRowNonZeroElementTable[ii].Node=node;

    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> * nodeTemp = node;
    int i0,j0,k0; //for test
    i0 = MatrixRowNonZeroElementTable[ii].i;
    j0 = MatrixRowNonZeroElementTable[ii].j;
    k0 = MatrixRowNonZeroElementTable[ii].k;

    if ((i0>=iMax) && (nodeTemp!=NULL)) {
      i0-=_BLOCK_CELLS_X_;
      nodeTemp=nodeTemp->GetNeibFace(1,0,0,PIC::Mesh::mesh);
    }
    else if (i0<0 && nodeTemp!=NULL) {
      i0+=_BLOCK_CELLS_X_;
      nodeTemp=nodeTemp->GetNeibFace(0,0,0,PIC::Mesh::mesh);
    }

    if ((j0>=jMax) && (nodeTemp!=NULL)) {
      j0-=_BLOCK_CELLS_Y_;
      nodeTemp=nodeTemp->GetNeibFace(3,0,0,PIC::Mesh::mesh);
    }
    else if (j0<0 && nodeTemp!=NULL) {
      j0+=_BLOCK_CELLS_Y_;
      nodeTemp=nodeTemp->GetNeibFace(2,0,0,PIC::Mesh::mesh);
    }

    if ((k0>=kMax) && (nodeTemp!=NULL)) {
      k0-=_BLOCK_CELLS_Z_;
      nodeTemp=nodeTemp->GetNeibFace(5,0,0,PIC::Mesh::mesh);
    }
    else if (k0<0 && nodeTemp!=NULL) {
      k0+=_BLOCK_CELLS_Z_;
      nodeTemp=nodeTemp->GetNeibFace(4,0,0,PIC::Mesh::mesh);
    }

    double xlocal[3];
    int indlocal[3]={MatrixRowNonZeroElementTable[ii].i,MatrixRowNonZeroElementTable[ii].j,MatrixRowNonZeroElementTable[ii].k};
    int indexG_local[3];
    bool isFixed=false;

    for (int idim=0; idim<3; idim++) {
      xlocal[idim]=MatrixRowNonZeroElementTable[ii].Node->xmin[idim]+indlocal[idim]*dx[idim];
    }
     
    pointLeft.push_back(ii);

    if (MatrixRowNonZeroElementTable[ii].Node==NULL){
      pointLeft.pop_back();
      continue;
    }
    else if (MatrixRowNonZeroElementTable[ii].Node->IsUsedInCalculationFlag==false) {
      pointLeft.pop_back();
      continue;
    }
    
    if (nodeTemp==NULL){
      pointLeft.pop_back();
      continue;
    }
    else if (nodeTemp->IsUsedInCalculationFlag==false){
      pointLeft.pop_back();
      continue;
    }

    if (isBoundaryCorner(xlocal,node)) pointLeft.pop_back();
  }

  for (int ii=0; ii<pointLeft.size();ii++){
    int copyFrom = pointLeft[ii];

    if (ii!=copyFrom){
      MatrixRowNonZeroElementTable[ii].i=MatrixRowNonZeroElementTable[copyFrom].i;
      MatrixRowNonZeroElementTable[ii].j=MatrixRowNonZeroElementTable[copyFrom].j;
      MatrixRowNonZeroElementTable[ii].k=MatrixRowNonZeroElementTable[copyFrom].k;

      MatrixRowNonZeroElementTable[ii].MatrixElementValue= MatrixRowNonZeroElementTable[copyFrom].MatrixElementValue;
      MatrixRowNonZeroElementTable[ii].iVar=MatrixRowNonZeroElementTable[copyFrom].iVar;

      MatrixRowNonZeroElementTable[ii].MatrixElementParameterTable[0]=MatrixRowNonZeroElementTable[copyFrom].MatrixElementParameterTable[0];
      MatrixRowNonZeroElementTable[ii].MatrixElementParameterTableLength=MatrixRowNonZeroElementTable[copyFrom].MatrixElementParameterTableLength;
      MatrixRowNonZeroElementTable[ii].MatrixElementSupportTableLength = MatrixRowNonZeroElementTable[copyFrom].MatrixElementSupportTableLength;
         
      MatrixRowNonZeroElementTable[ii].MatrixElementSupportTable[0]=MatrixRowNonZeroElementTable[copyFrom].MatrixElementSupportTable[0];
    }

  }
  
  NonZeroElementsFound=pointLeft.size();
  
  //NonZeroElementsFound=81;

  //  for (int iVarIndex=0; iVarIndex<3; iVarIndex++){

    // fill first 27 elements
  for (int ii=0;ii<3;ii++){
    for (int jj=0;jj<3;jj++){
      for (int kk=0;kk<3;kk++){
        int iElement = ii+jj*3+kk*9;
        int iNode = i+indexAddition[ii];
        int jNode = j+indexAddition[jj];
        int kNode = k+indexAddition[kk];

      {
        const int di = indexAddition[ii];
        const int dj = indexAddition[jj];
        const int dk = indexAddition[kk];
        //indexer[iVar].check_and_set(iElement,di, dj, dk, 0);
	iElement=indexer[iVar].get_idx(di, dj, dk, 0);
      }

        RhsSupportTable_CornerNodes[iElement].Coefficient= 0.0;
        RhsSupportTable_CornerNodes[iElement].AssociatedDataPointer=node->block->GetCornerNode(_getCornerNodeLocalNumber(iNode,jNode,kNode))->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset;

        // Store a semantic RHS support entry: Corner/E/component=0 (Ex).
        auto *cornerNodePtr = node->block->GetCornerNode(_getCornerNodeLocalNumber(iNode,jNode,kNode));
        char *pntTemp = (cornerNodePtr!=NULL) ? cornerNodePtr->GetAssociatedDataBufferPointer() : NULL;

        RhsSupportTable_CornerNodes[iElement].SetCornerE(0.0, (pntTemp!=NULL) ? cornerNodePtr : NULL, 0);

      }
    }
  }


  for (int iVarIndex=1; iVarIndex<3; iVarIndex++){
    // fill next 54 elements
    for (int ii=0;ii<3;ii++){
      for (int jj=0;jj<3;jj++){
        for (int kk=0;kk<3;kk++){

          int iElement = iVarIndex*27+ii+jj*3+kk*9;
          int jOldElement = ii+jj*3+kk*9;

      {
        const int di = indexOffset[ii];
        const int dj = indexOffset[jj];
        const int dk = indexOffset[kk];
        //indexer[iVar].check_and_set(iElement,di, dj, dk, iVarIndex);
	iElement=indexer[iVar].get_idx(di, dj, dk, iVarIndex);

        //indexer[iVar].check_and_set(jOldElement,di, dj, dk, 0);
	jOldElement=indexer[iVar].get_idx(di, dj, dk, 0);

      }


          RhsSupportTable_CornerNodes[iElement].Coefficient= 0.0;
          RhsSupportTable_CornerNodes[iElement].AssociatedDataPointer=RhsSupportTable_CornerNodes[jOldElement].AssociatedDataPointer;

	            RhsSupportTable_CornerNodes[iElement].SetCornerE(0.0, RhsSupportTable_CornerNodes[jOldElement].corner, (unsigned char)iVarIndex);
        }
      }
    }
  }


  //  bool isTest=false;
  //if (fabs(x[0]-1.0)<0.1 && fabs(x[1]-1.0)<0.1 && fabs(x[2]-1.0)<0.1)
  //  isTest = true;
if (_PIC_STENCIL_NUMBER_==375) {  
  for (int iVarIndex=0; iVarIndex<1; iVarIndex++){
    int cntTemp = 0;

    for (int kk=0; kk<5; kk++){
      for (int jj=0; jj<5; jj++){
        for (int ii=0; ii<5; ii++){

          if (ii<3 && jj<3 && kk<3) continue;

          int iNode = i+indexOffset[ii];
          int jNode = j+indexOffset[jj];
          int kNode = k+indexOffset[kk];

          int iElement = 81+iVarIndex*98+cntTemp;
          cntTemp++;


      {
        const int di = indexOffset[ii];
        const int dj = indexOffset[jj];
        const int dk = indexOffset[kk];
        //indexer[iVar].check_and_set(iElement,di, dj, dk, iVarIndex);
	iElement=indexer[iVar].get_idx(di, dj, dk, iVarIndex);
      }


          RhsSupportTable_CornerNodes[iElement].Coefficient= 0.0;




          char * pntTemp = node->block->GetCornerNode(_getCornerNodeLocalNumber(iNode,jNode,kNode))->GetAssociatedDataBufferPointer();


	  auto *cornerNodePtr = node->block->GetCornerNode(_getCornerNodeLocalNumber(iNode,jNode,kNode));
	            RhsSupportTable_CornerNodes[iElement].SetCornerE(0.0, (pntTemp!=NULL) ? cornerNodePtr : NULL, 0);

          if (pntTemp) {
            RhsSupportTable_CornerNodes[iElement].AssociatedDataPointer=pntTemp+PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset;
          }
          else{
            RhsSupportTable_CornerNodes[iElement].AssociatedDataPointer = NULL;
          }
        }
      }
    }
  }


  for (int iVarIndex=1; iVarIndex<3; iVarIndex++){
    int cntTemp = 0;
    for (int kk=0; kk<5; kk++){
      for (int jj=0; jj<5; jj++){
        for (int ii=0; ii<5; ii++){
          if (ii<3 && jj<3 && kk<3) continue;

          int iElement = 81+iVarIndex*98+cntTemp;
          int iElementOld = 81+cntTemp;


      {
        const int di = indexOffset[ii];
        const int dj = indexOffset[jj];
        const int dk = indexOffset[kk];
        //indexer[iVar].check_and_set(iElement,di, dj, dk, iVarIndex);
	iElement=indexer[iVar].get_idx(di, dj, dk, iVarIndex);

        //indexer[iVar].check_and_set(iElementOld,di, dj, dk, 0);
	iElementOld=indexer[iVar].get_idx(di, dj, dk, 0);

      }


          cntTemp++;

          RhsSupportTable_CornerNodes[iElement].Coefficient= 0.0;
          RhsSupportTable_CornerNodes[iElement].AssociatedDataPointer=RhsSupportTable_CornerNodes[iElementOld].AssociatedDataPointer;


	              RhsSupportTable_CornerNodes[iElement].SetCornerE(0.0, RhsSupportTable_CornerNodes[iElementOld].corner, (unsigned char)iVarIndex);

        }
      }
    }
  }
}
    

  //laplacian
  for (int idim=0;idim<3;idim++) {
    cStencil::cStencilData *st=LaplacianStencil+idim;

    for (int it=0;it<st->Length;it++) {
      int ii=reversed_indexAddition[st->Data[it].i+1];
      int jj=reversed_indexAddition[st->Data[it].j+1];
      int kk=reversed_indexAddition[st->Data[it].k+1];

      int index=ii+jj*3+kk*9;
      int iElement = index + iVar*27;

      {
        const int di = st->Data[it].i;
        const int dj = st->Data[it].j;
        const int dk = st->Data[it].k;
        //indexer[iVar].check_and_set(iElement,di, dj, dk, iVar);
	iElement=indexer[iVar].get_idx(di, dj, dk, iVar);
      }

      //plus laplacian
      RhsSupportTable_CornerNodes[iElement].Coefficient+=st->Data[it].a*coeffSqr[idim];

      RhsSupportTable_CornerNodes[iElement].CoefficientNEW+=st->Data[it].a*coeffSqr[idim];
    }
  }


  //minus graddiv E
  for (int iVarIndex=0;iVarIndex<3;iVarIndex++){
    cStencil::cStencilData *st=&GradDivStencil[iVar][iVarIndex];


    for (int it=0;it<st->Length;it++) {
      int ii=reversed_indexAddition[st->Data[it].i+1];
      int jj=reversed_indexAddition[st->Data[it].j+1];
      int kk=reversed_indexAddition[st->Data[it].k+1];

      int nodeIndex=ii+jj*3+kk*9;
      int iElement = nodeIndex + iVarIndex*27;

      {
        const int di = st->Data[it].i;
        const int dj = st->Data[it].j;
        const int dk = st->Data[it].k;
        //indexer[iVar].check_and_set(iElement,di, dj, dk, iVarIndex);
	iElement=indexer[iVar].get_idx(di, dj, dk, iVarIndex);
      }

      RhsSupportTable_CornerNodes[iElement].Coefficient -=(1-corrCoeff)*st->Data[it].a*coeff[iVar]*coeff[iVarIndex];

      RhsSupportTable_CornerNodes[iElement].CoefficientNEW -=(1-corrCoeff)*st->Data[it].a*coeff[iVar]*coeff[iVarIndex];
    }
  }

if (_PIC_STENCIL_NUMBER_==375) {
  for (int iVarIndex=0;iVarIndex<3;iVarIndex++){
    cStencil::cStencilData *st375=&GradDivStencil375[iVar][iVarIndex];

    for (int it=0;it<st375->Length;it++) {

      if ((st375->Data[it].i<-1)||(st375->Data[it].i>1) || (st375->Data[it].j<-1)||(st375->Data[it].j>1) ||(st375->Data[it].k<-1)||(st375->Data[it].k>1) ) continue;

      int ii=reversed_indexAddition[st375->Data[it].i+1];
      int jj=reversed_indexAddition[st375->Data[it].j+1];
      int kk=reversed_indexAddition[st375->Data[it].k+1];

      int nodeIndex=ii+jj*3+kk*9;
      int iElement = nodeIndex + iVarIndex*27;

      {
        const int di = st375->Data[it].i;
        const int dj = st375->Data[it].j;
        const int dk = st375->Data[it].k;
        //indexer[iVar].check_and_set(iElement,di, dj, dk, iVarIndex);
	iElement=indexer[iVar].get_idx(di, dj, dk, iVarIndex);
      }


      RhsSupportTable_CornerNodes[iElement].Coefficient -=corrCoeff*st375->Data[it].a*coeff[iVar]*coeff[iVarIndex];
      RhsSupportTable_CornerNodes[iElement].CoefficientNEW -=corrCoeff*st375->Data[it].a*coeff[iVar]*coeff[iVarIndex];
    }
  }


  for (int iVarIndex=0; iVarIndex<3; iVarIndex++){
    cStencil::cStencilData *st375=&GradDivStencil375[iVar][iVarIndex];

    for (int it=0;it<st375->Length;it++) {
      int ii=reversed_indexOffset[st375->Data[it].i+2];
      int jj=reversed_indexOffset[st375->Data[it].j+2];
      int kk=reversed_indexOffset[st375->Data[it].k+2];

      if (ii<3 && jj<3 && kk<3) continue;

      int iElement = 81+iVarIndex*98+OrderingOffsetTable[ii][jj][kk];
      int nodeIndex = 27+OrderingOffsetTable[ii][jj][kk];

      {
        const int di = st375->Data[it].i;
        const int dj = st375->Data[it].j;
        const int dk = st375->Data[it].k;
        //indexer[iVar].check_and_set(iElement,di, dj, dk, iVarIndex);
	iElement=indexer[iVar].get_idx(di, dj, dk, iVarIndex);
      }


      RhsSupportTable_CornerNodes[iElement].Coefficient -=corrCoeff*st375->Data[it].a*coeff[iVar]*coeff[iVarIndex];
      RhsSupportTable_CornerNodes[iElement].CoefficientNEW -=corrCoeff*st375->Data[it].a*coeff[iVar]*coeff[iVarIndex];
    }
  }

}

int idx=indexer[iVar].get_index_count(); //_PIC_STENCIL_NUMBER_  
int idx2;
      {
        const int di = 0;
        const int dj = 0;
        const int dk = 0;
        //indexer[iVar].check_and_set(27*iVar,di, dj, dk, iVar);
	idx2=indexer[iVar].get_idx(di, dj, dk, iVar);
      }

  RhsSupportTable_CornerNodes[idx].AssociatedDataPointer=RhsSupportTable_CornerNodes[idx2].AssociatedDataPointer;
  RhsSupportTable_CornerNodes[idx].Coefficient=-4*Pi*dtTotal*theta;

  RhsSupportTable_CornerNodes[idx].SetCornerJ(-4*Pi*dtTotal*theta, RhsSupportTable_CornerNodes[idx2].corner, (unsigned char)iVar);
 
  RhsSupportLength_CornerNodes=idx+1;

  //Ex^n,Ey^n,Ez^n
  rhs=0.0;

  int indexAdditionB[2] = {-1,0};
  int iElement = 0;
  
  double curlB = 0.0;
    //Ex  rhs+= d Bz/dy - d By/dz
  if (iVar==0){
          
      for (int ii=0;ii<2;ii++){
        for (int jj=0;jj<2;jj++){

        PIC::Mesh::cDataCenterNode* center = node->block->GetCenterNode(
            _getCenterNodeLocalNumber(i+indexAdditionB[ii], j, k+indexAdditionB[jj]));
        
        RhsSupportTable_CenterNodes[iElement].SetCenterB(
            coeff4[1],    // +c*dt/(2*dy)
            center,
            0             // component 2 = Bz
        );

          RhsSupportTable_CenterNodes[iElement].Coefficient=coeff4[1]; //c(dt)/dy
          RhsSupportTable_CenterNodes[iElement].AssociatedDataPointer=node->block->GetCenterNode(_getCenterNodeLocalNumber(i+indexAdditionB[ii],j,k+indexAdditionB[jj]))->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset;
          //  rhs+=((double*)(RhsSupportTable_CenterNodes[iElement].AssociatedDataPointer+CurrentCenterNodeOffset))[BzOffsetIndex]*RhsSupportTable_CenterNodes[iElement].Coefficient;
          iElement++;

        }
      }

      for (int ii=0;ii<2;ii++){
        for (int jj=0;jj<2;jj++){

        PIC::Mesh::cDataCenterNode* center = node->block->GetCenterNode(
            _getCenterNodeLocalNumber(i+indexAdditionB[ii], j-1, k+indexAdditionB[jj]));
        
        RhsSupportTable_CenterNodes[iElement].SetCenterB(
            -coeff4[1],   // -c*dt/(2*dy)
            center,
            0             // component 2 = Bz
        );

          RhsSupportTable_CenterNodes[iElement].Coefficient=-coeff4[1]; //-c(dt)/dy
          RhsSupportTable_CenterNodes[iElement].AssociatedDataPointer=node->block->GetCenterNode(_getCenterNodeLocalNumber(i+indexAdditionB[ii],j-1,k+indexAdditionB[jj]))->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset;
          // rhs+=((double*)(RhsSupportTable_CenterNodes[iElement].AssociatedDataPointer+CurrentCenterNodeOffset))[BzOffsetIndex]*RhsSupportTable_CenterNodes[iElement].Coefficient;
          iElement++;
        }
      }
      
      for (int ii=0;ii<2;ii++){
        for (int jj=0;jj<2;jj++){
        PIC::Mesh::cDataCenterNode* center = node->block->GetCenterNode(
            _getCenterNodeLocalNumber(i+indexAdditionB[ii], j+indexAdditionB[jj], k));
        
        RhsSupportTable_CenterNodes[iElement].SetCenterB(
            -coeff4[2],   // -c*dt/(2*dz)
            center,
            0             // component 1 = By
        );

          RhsSupportTable_CenterNodes[iElement].Coefficient=-coeff4[2]; //-c(dt)/dz
          RhsSupportTable_CenterNodes[iElement].AssociatedDataPointer=node->block->GetCenterNode(_getCenterNodeLocalNumber(i+indexAdditionB[ii],j+indexAdditionB[jj],k))->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset;
          //  rhs+=((double*)(RhsSupportTable_CenterNodes[iElement].AssociatedDataPointer+CurrentCenterNodeOffset))[ByOffsetIndex]*RhsSupportTable_CenterNodes[iElement].Coefficient;
          iElement++;
        }
      }
    
      for (int ii=0;ii<2;ii++){
        for (int jj=0;jj<2;jj++){

        PIC::Mesh::cDataCenterNode* center = node->block->GetCenterNode(
            _getCenterNodeLocalNumber(i+indexAdditionB[ii], j+indexAdditionB[jj], k-1));
        
        RhsSupportTable_CenterNodes[iElement].SetCenterB(
            coeff4[2],    // +c*dt/(2*dz)
            center,
            0             // component 1 = By
        );

          RhsSupportTable_CenterNodes[iElement].Coefficient=coeff4[2]; //c(dt)/dz
          RhsSupportTable_CenterNodes[iElement].AssociatedDataPointer=node->block->GetCenterNode(_getCenterNodeLocalNumber(i+indexAdditionB[ii],j+indexAdditionB[jj],k-1))->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset;
          // rhs+=((double*)(RhsSupportTable_CenterNodes[iElement].AssociatedDataPointer+CurrentCenterNodeOffset))[ByOffsetIndex]*RhsSupportTable_CenterNodes[iElement].Coefficient;
          iElement++;
        }
      }
    
    }

     //Ey  rhs+= d Bx/dz - d Bz/dx
    if (iVar==1){
      for (int ii=0;ii<2;ii++){
        for (int jj=0;jj<2;jj++){
        PIC::Mesh::cDataCenterNode* center = node->block->GetCenterNode(
            _getCenterNodeLocalNumber(i+indexAdditionB[ii], j+indexAdditionB[jj], k));
        
        RhsSupportTable_CenterNodes[iElement].SetCenterB(
            coeff4[2],    // +c*dt/(2*dz)
            center,
            0             // component 0 = Bx
        );

          RhsSupportTable_CenterNodes[iElement].Coefficient=coeff4[2]; //c(dt)/dz
          RhsSupportTable_CenterNodes[iElement].AssociatedDataPointer=node->block->GetCenterNode(_getCenterNodeLocalNumber(i+indexAdditionB[ii],j+indexAdditionB[jj],k))->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset;
          // rhs+=((double*)(RhsSupportTable_CenterNodes[iElement].AssociatedDataPointer+CurrentCenterNodeOffset))[BxOffsetIndex]*RhsSupportTable_CenterNodes[iElement].Coefficient;
          // curlB+=((double*)(RhsSupportTable_CenterNodes[iElement].AssociatedDataPointer+CurrentCenterNodeOffset))[BxOffsetIndex]*RhsSupportTable_CenterNodes[iElement].Coefficient;
          iElement++;

        }
      }

      for (int ii=0;ii<2;ii++){
        for (int jj=0;jj<2;jj++){

        PIC::Mesh::cDataCenterNode* center = node->block->GetCenterNode(
            _getCenterNodeLocalNumber(i+indexAdditionB[ii], j+indexAdditionB[jj], k-1));
        
        RhsSupportTable_CenterNodes[iElement].SetCenterB(
            -coeff4[2],   // -c*dt/(2*dz)
            center,
            0             // component 0 = Bx
        );

          RhsSupportTable_CenterNodes[iElement].Coefficient=-coeff4[2]; //-c(dt)/dz
          RhsSupportTable_CenterNodes[iElement].AssociatedDataPointer=node->block->GetCenterNode(_getCenterNodeLocalNumber(i+indexAdditionB[ii],j+indexAdditionB[jj],k-1))->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset;
          //rhs+=((double*)(RhsSupportTable_CenterNodes[iElement].AssociatedDataPointer+CurrentCenterNodeOffset))[BxOffsetIndex]*RhsSupportTable_CenterNodes[iElement].Coefficient;
          //curlB+=((double*)(RhsSupportTable_CenterNodes[iElement].AssociatedDataPointer+CurrentCenterNodeOffset))[BxOffsetIndex]*RhsSupportTable_CenterNodes[iElement].Coefficient;
          iElement++;
        }
      }
      
      for (int ii=0;ii<2;ii++){
        for (int jj=0;jj<2;jj++){

        PIC::Mesh::cDataCenterNode* center = node->block->GetCenterNode(
            _getCenterNodeLocalNumber(i, j+indexAdditionB[jj], k+indexAdditionB[ii]));
        
        RhsSupportTable_CenterNodes[iElement].SetCenterB(
            -coeff4[0],   // -c*dt/(2*dx)
            center,
            0             // component 2 = Bz
        );          

          RhsSupportTable_CenterNodes[iElement].Coefficient=-coeff4[0]; //-c(dt)/dx
          RhsSupportTable_CenterNodes[iElement].AssociatedDataPointer=node->block->GetCenterNode(_getCenterNodeLocalNumber(i,j+indexAdditionB[jj],k+indexAdditionB[ii]))->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset;
          //rhs+=((double*)(RhsSupportTable_CenterNodes[iElement].AssociatedDataPointer+CurrentCenterNodeOffset))[BzOffsetIndex]*RhsSupportTable_CenterNodes[iElement].Coefficient;
          //curlB+=((double*)(RhsSupportTable_CenterNodes[iElement].AssociatedDataPointer+CurrentCenterNodeOffset))[BzOffsetIndex]*RhsSupportTable_CenterNodes[iElement].Coefficient;
          iElement++;
        }
      }
    
      for (int ii=0;ii<2;ii++){
        for (int jj=0;jj<2;jj++){

        PIC::Mesh::cDataCenterNode* center = node->block->GetCenterNode(
            _getCenterNodeLocalNumber(i-1, j+indexAdditionB[jj], k+indexAdditionB[ii]));
        
        RhsSupportTable_CenterNodes[iElement].SetCenterB(
            coeff4[0],    // +c*dt/(2*dx)
            center,
            0             // component 2 = Bz
        );

          RhsSupportTable_CenterNodes[iElement].Coefficient=coeff4[0]; //c(dt)/dx
          RhsSupportTable_CenterNodes[iElement].AssociatedDataPointer=node->block->GetCenterNode(_getCenterNodeLocalNumber(i-1,j+indexAdditionB[jj],k+indexAdditionB[ii]))->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset;
          // rhs+=((double*)(RhsSupportTable_CenterNodes[iElement].AssociatedDataPointer+CurrentCenterNodeOffset))[BzOffsetIndex]*RhsSupportTable_CenterNodes[iElement].Coefficient;
          //curlB+=((double*)(RhsSupportTable_CenterNodes[iElement].AssociatedDataPointer+CurrentCenterNodeOffset))[BzOffsetIndex]*RhsSupportTable_CenterNodes[iElement].Coefficient;
          iElement++;
        }
      }
    }
    
    //Ez  rhs+= d By/dx - d Bx/dy
    if (iVar==2) {
      for (int ii=0;ii<2;ii++){
        for (int jj=0;jj<2;jj++){
        PIC::Mesh::cDataCenterNode* center = node->block->GetCenterNode(
            _getCenterNodeLocalNumber(i, j+indexAdditionB[jj], k+indexAdditionB[ii]));
        
        RhsSupportTable_CenterNodes[iElement].SetCenterB(
            coeff4[0],    // +c*dt/(2*dx)
            center,
            0             // component 1 = By
        );


          RhsSupportTable_CenterNodes[iElement].Coefficient=coeff4[0]; //c(dt)/dx
          RhsSupportTable_CenterNodes[iElement].AssociatedDataPointer=node->block->GetCenterNode(_getCenterNodeLocalNumber(i,j+indexAdditionB[jj],k+indexAdditionB[ii]))->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset;
          //rhs+=((double*)(RhsSupportTable_CenterNodes[iElement].AssociatedDataPointer+CurrentCenterNodeOffset))[ByOffsetIndex]*RhsSupportTable_CenterNodes[ii].Coefficient;
          //curlB+=((double*)(RhsSupportTable_CenterNodes[iElement].AssociatedDataPointer+CurrentCenterNodeOffset))[ByOffsetIndex]*RhsSupportTable_CenterNodes[ii].Coefficient;
          iElement++;
        }
      }

      for (int ii=0;ii<2;ii++){
        for (int jj=0;jj<2;jj++){

        PIC::Mesh::cDataCenterNode* center = node->block->GetCenterNode(
            _getCenterNodeLocalNumber(i-1, j+indexAdditionB[jj], k+indexAdditionB[ii]));
        
        RhsSupportTable_CenterNodes[iElement].SetCenterB(
            -coeff4[0],   // -c*dt/(2*dx)
            center,
            0             // component 1 = By
        );

          RhsSupportTable_CenterNodes[iElement].Coefficient=-coeff4[0]; //-c(dt)/dx
          RhsSupportTable_CenterNodes[iElement].AssociatedDataPointer=node->block->GetCenterNode(_getCenterNodeLocalNumber(i-1,j+indexAdditionB[jj],k+indexAdditionB[ii]))->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset;
          //rhs+=((double*)(RhsSupportTable_CenterNodes[iElement].AssociatedDataPointer+CurrentCenterNodeOffset))[ByOffsetIndex]*RhsSupportTable_CenterNodes[iElement].Coefficient;
          //curlB+=((double*)(RhsSupportTable_CenterNodes[iElement].AssociatedDataPointer+CurrentCenterNodeOffset))[ByOffsetIndex]*RhsSupportTable_CenterNodes[iElement].Coefficient;
          iElement++;
        }
      }
      
      for (int ii=0;ii<2;ii++){
        for (int jj=0;jj<2;jj++){

        PIC::Mesh::cDataCenterNode* center = node->block->GetCenterNode(
            _getCenterNodeLocalNumber(i+indexAdditionB[jj], j, k+indexAdditionB[ii]));
        
        RhsSupportTable_CenterNodes[iElement].SetCenterB(
            -coeff4[1],   // -c*dt/(2*dy)
            center,
            0             // component 0 = Bx
        );

          RhsSupportTable_CenterNodes[iElement].Coefficient=-coeff4[1]; //-c(dt)/dy
          RhsSupportTable_CenterNodes[iElement].AssociatedDataPointer=node->block->GetCenterNode(_getCenterNodeLocalNumber(i+indexAdditionB[jj],j,k+indexAdditionB[ii]))->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset;
          //rhs+=((double*)(RhsSupportTable_CenterNodes[iElement].AssociatedDataPointer+CurrentCenterNodeOffset))[BxOffsetIndex]*RhsSupportTable_CenterNodes[iElement].Coefficient;
          //curlB+=((double*)(RhsSupportTable_CenterNodes[iElement].AssociatedDataPointer+CurrentCenterNodeOffset))[BxOffsetIndex]*RhsSupportTable_CenterNodes[iElement].Coefficient;
          iElement++;
        }
      }
    
      for (int ii=0;ii<2;ii++){
        for (int jj=0;jj<2;jj++){

        PIC::Mesh::cDataCenterNode* center = node->block->GetCenterNode(
            _getCenterNodeLocalNumber(i+indexAdditionB[jj], j-1, k+indexAdditionB[ii]));
        
        RhsSupportTable_CenterNodes[iElement].SetCenterB(
            coeff4[1],    // +c*dt/(2*dy)
            center,
            0             // component 0 = Bx
        );
                  
          RhsSupportTable_CenterNodes[iElement].Coefficient=coeff4[1]; //c(dt)/dy
          RhsSupportTable_CenterNodes[iElement].AssociatedDataPointer=node->block->GetCenterNode(_getCenterNodeLocalNumber(i+indexAdditionB[jj],j-1,k+indexAdditionB[ii]))->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset;
          //rhs+=((double*)(RhsSupportTable_CenterNodes[iElement].AssociatedDataPointer+CurrentCenterNodeOffset))[BxOffsetIndex]*RhsSupportTable_CenterNodes[iElement].Coefficient;
          //curlB+=((double*)(RhsSupportTable_CenterNodes[iElement].AssociatedDataPointer+CurrentCenterNodeOffset))[BxOffsetIndex]*RhsSupportTable_CenterNodes[iElement].Coefficient;
          iElement++;
        }
      }

      //double analytic = -1000*3.14159265/2*cos((x[0]+1)*3.14159265/2)*0.2;
      //printf("Ez,curlB:%f,analytic:%f\n", curlB, analytic);
      //rhs+=curlB;
    }
   
    RhsSupportLength_CenterNodes = iElement;     

if (false) { 
  PostAssembleSelfCheck(i, j, k, iVar,
                        MatrixRowNonZeroElementTable, NonZeroElementsFound, rhs,
                        RhsSupportTable_CornerNodes, RhsSupportLength_CornerNodes,
                        RhsSupportTable_CenterNodes, RhsSupportLength_CenterNodes,
                        node);
}

//GetStencil(i,j,k,iVar,MatrixRowNonZeroElementTable,NonZeroElementsFound,rhs,RhsSupportTable_CornerNodes,RhsSupportLength_CornerNodes,RhsSupportTable_CenterNodes,RhsSupportLength_CenterNodes,node);
}


