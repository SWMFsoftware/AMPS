//$Id$
//interface to SWMF's GMRES solver

/*
 * LinearSystemCornerNode.h
 *
 *  Created on: Nov 20, 2017
 *      Author: vtenishe
 */

#include "linear_solver_wrapper_c.h"
#include <functional>
#include <iostream>

#if _AVX_INSTRUCTIONS_USAGE_MODE_ == _AVX_INSTRUCTIONS_USAGE_MODE__ON_
#include <immintrin.h>
#endif

#ifndef _LINEARSYSTEMCORNERNODE_H_
#define _LINEARSYSTEMCORNERNODE_H_

#ifndef _PIC_MATMUL_MPI_MULTITHREAD_
#define _PIC_MATMUL_MPI_MULTITHREAD_ _PIC_MODE_OFF_
#endif



#if _AVX_INSTRUCTIONS_USAGE_MODE_ == _AVX_INSTRUCTIONS_USAGE_MODE__512_
  #ifndef _AVX_MATMUL_ 
  #define _AVX_MATMUL_ _ON_
  #endif
#elif _AVX_INSTRUCTIONS_USAGE_MODE_ == _AVX_INSTRUCTIONS_USAGE_MODE__256_
  #ifndef _AVX_MATMUL_ 
  #define _AVX_MATMUL_ _ON_
  #endif
#else 
  #ifdef _AVX_MATMUL_
  #undef _AVX_MATMUL_
  #endif

  #define  _AVX_MATMUL_  _OFF_
#endif

class cLinearSystemCornerNodeDataRequestListElement {
public:
  int CornerNodeID;
  cAMRnodeID NodeID;
};

class cRebuildMatrix {
public:
  virtual void RebuildMatrix() = 0;
};

template <class cCornerNode, int NodeUnknownVariableVectorLength,int MaxStencilLength,
int MaxRhsSupportLength_CornerNodes,int MaxRhsSupportLength_CenterNodes,
int MaxMatrixElementParameterTableLength,int MaxMatrixElementSupportTableLength>
class cLinearSystemCornerNode : public cRebuildMatrix {
public:
  double *SubdomainPartialRHS,*SubdomainPartialUnknownsVector;
  int SubdomainPartialUnknownsVectorLength;

  //counter of the matrix modifications 
  int nMartixModifications;

  class cStencilElementData {
  public: 
    int UnknownVectorIndex,iVar,Thread;
    double MatrixElementValue;
    double *rhs;
  
    cStencilElementData() {
      UnknownVectorIndex=-1,iVar=-1,Thread=-1;
      MatrixElementValue=0.0;
      rhs=NULL;
    }


    bool operator==(const cStencilElementData& other) const {
      // Compare all integer members
      if (UnknownVectorIndex != other.UnknownVectorIndex) {
        return false;
      }

      if (iVar != other.iVar) {
        return false;
      }

      if (Thread != other.Thread) {
        return false;
      }

      // Compare double value for exact equality
      if (MatrixElementValue != other.MatrixElementValue) {
        return false;
      }

      // Compare pointer identity for rhs
      if (rhs != other.rhs) {
        return false;
      }

      return true;
    }

   // Inequality operator (optional but recommended)
   bool operator!=(const cStencilElementData& other) const {
     return !(*this == other);
   }

    // In cStencilElementData (public:), include <cstdio> somewhere in the TU.
    static bool Compare(const cStencilElementData& a,
                    const cStencilElementData& b,
                    bool break_flag = false) {
      bool equal = true;

      // Print header lazily (only once if something differs)
      auto print_header = [&]() {
        static bool printed = false;
        if (!printed) {
          if (break_flag) std::printf("cStencilElementData mismatch:\n");
          printed = true;
        }
      };

      if (a.UnknownVectorIndex != b.UnknownVectorIndex) {
       equal = false; print_header();
       if (break_flag) std::printf("  [UnknownVectorIndex] a=%d b=%d\n",
                a.UnknownVectorIndex, b.UnknownVectorIndex);
      }

      if (a.iVar != b.iVar) {
        equal = false; print_header();
        if (break_flag) std::printf("  [iVar] a=%d b=%d\n", a.iVar, b.iVar);
      }

      if (a.Thread != b.Thread) {
        equal = false; print_header();
        if (break_flag) std::printf("  [Thread] a=%d b=%d\n", a.Thread, b.Thread);
      }

      if (a.MatrixElementValue != b.MatrixElementValue) {
        equal = false; print_header();
        if (break_flag) std::printf("  [MatrixElementValue] a=%.17e b=%.17e\n",
                a.MatrixElementValue, b.MatrixElementValue);
      }

      if (a.rhs != b.rhs) {
        equal = false; print_header();
        if (break_flag) std::printf("  [rhs ptr] a=%p b=%p\n",
                (const void*)a.rhs, (const void*)b.rhs);
      }

      if (!equal && break_flag) {
        // Build a short message for the project-specific exit() hook.
        char msg[256];
        std::snprintf(msg, sizeof(msg),
                  "cStencilElementData objects differ (see printf log above)");
        exit(__LINE__, __FILE__, msg);
      }

      return equal;
    }

  };

  class cStencilElement {
  public:
    cCornerNode* CornerNode;
    int CornerNodeID;
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;

    //the user-model parameters neede for recalculating of the matrix coefficients
    double MatrixElementParameterTable[MaxMatrixElementParameterTableLength];
    int MatrixElementParameterTableLength;

    //the user-defined support used for recalculating the matrix coefficients
    void *MatrixElementSupportTable[MaxMatrixElementSupportTableLength];
    int MatrixElementSupportTableLength;

    cStencilElement()  {
      node=NULL,CornerNode=NULL;
      CornerNodeID=-1;

      MatrixElementParameterTableLength=0,MatrixElementSupportTableLength=0;

      for (int i=0;i<MaxMatrixElementParameterTableLength;i++) MatrixElementParameterTable[i]=0.0;
      for (int i=0;i<MaxMatrixElementSupportTableLength;i++) MatrixElementSupportTable[i]=NULL;
    }

    // EQUALITY OPERATORS (already provided earlier)
    bool operator==(const cStencilElement& other) const {
      if (CornerNode != other.CornerNode || node != other.node) {
        return false;
      }

      if (CornerNodeID != other.CornerNodeID ||
        MatrixElementParameterTableLength != other.MatrixElementParameterTableLength ||
        MatrixElementSupportTableLength != other.MatrixElementSupportTableLength) {
        return false;
      }

      for (int i = 0; i < MatrixElementParameterTableLength; ++i) {
        if (MatrixElementParameterTable[i] != other.MatrixElementParameterTable[i]) {
            return false;
        }
      }

      for (int i = 0; i < MatrixElementSupportTableLength; ++i) {
        if (MatrixElementSupportTable[i] != other.MatrixElementSupportTable[i]) {
            return false;
        }
      }

      return true;
   }

   bool operator!=(const cStencilElement& other) const {
    return !(*this == other);
   }

    static bool Compare(const cStencilElement& a,
                    const cStencilElement& b,
                    bool break_flag = false) {
      bool equal = true;
      bool printed_header = false;

      auto hdr = [&](){
        if (!printed_header) { std::printf("cStencilElement mismatch:\n"); printed_header = true; }
      };

      // --- Strict pointer identity checks (must match) ---
      if (a.CornerNode != b.CornerNode) {
        equal = false; hdr();
        if (break_flag) std::printf("  [CornerNode ptr] a=%p b=%p\n",
                (const void*)a.CornerNode, (const void*)b.CornerNode);
      }
      if (a.node != b.node) {
        equal = false; hdr();
        if (break_flag) std::printf("  [node ptr] a=%p b=%p\n",
                (const void*)a.node, (const void*)b.node);
      }
      if (a.MatrixElementParameterTable != b.MatrixElementParameterTable) {
        equal = false; hdr();
        if (break_flag) std::printf("  [ParamTable ptr] a=%p b=%p\n",
                (const void*)a.MatrixElementParameterTable,
                (const void*)b.MatrixElementParameterTable);
      }
      if (a.MatrixElementSupportTable != b.MatrixElementSupportTable) {
        equal = false; hdr();
        if (break_flag) std::printf("  [SupportTable ptr] a=%p b=%p\n",
                (const void*)a.MatrixElementSupportTable,
                (const void*)b.MatrixElementSupportTable);
      }

      // IDs/lengths must match too.
      if (a.CornerNodeID != b.CornerNodeID) {
        equal = false; hdr();
        if (break_flag) std::printf("  [CornerNodeID] a=%d b=%d\n", a.CornerNodeID, b.CornerNodeID);
      }
      if (a.MatrixElementParameterTableLength != b.MatrixElementParameterTableLength) {
        equal = false; hdr();
        if (break_flag) std::printf("  [ParamLen] a=%d b=%d\n",
                a.MatrixElementParameterTableLength, b.MatrixElementParameterTableLength);
      }
      if (a.MatrixElementSupportTableLength != b.MatrixElementSupportTableLength) {
        equal = false; hdr();
        if (break_flag) std::printf("  [SupportLen] a=%d b=%d\n",
                a.MatrixElementSupportTableLength, b.MatrixElementSupportTableLength);
      }

      // Optional: if you also want to flag per-entry differences (useful when pointers match
      // but contents might have changed in-place). This does NOT relax the pointer rule above.
      {
        const int len_min = (a.MatrixElementParameterTableLength < b.MatrixElementParameterTableLength)
                        ? a.MatrixElementParameterTableLength : b.MatrixElementParameterTableLength;
        for (int i = 0; i < len_min; ++i) {
          const double av = a.MatrixElementParameterTable ? a.MatrixElementParameterTable[i] : 0.0;
          const double bv = b.MatrixElementParameterTable ? b.MatrixElementParameterTable[i] : 0.0;
          if (av != bv) {
            equal = false; hdr();
            if (break_flag) std::printf("  [Param[%d]] a=%.17e b=%.17e\n", i, av, bv);
          }
        }
      }
      {
        const int len_min = (a.MatrixElementSupportTableLength < b.MatrixElementSupportTableLength)
                        ? a.MatrixElementSupportTableLength : b.MatrixElementSupportTableLength;
        for (int i = 0; i < len_min; ++i) {
          const void* ap = a.MatrixElementSupportTable ? a.MatrixElementSupportTable[i] : nullptr;
          const void* bp = b.MatrixElementSupportTable ? b.MatrixElementSupportTable[i] : nullptr;

          if (ap != bp) {
            equal = false; hdr();
            if (break_flag) std::printf("  [Support[%d] ptr] a=%p b=%p\n", i, ap, bp);
          }
        }
      }

      if (!equal && break_flag) {
        char msg[256];
        std::snprintf(msg, sizeof(msg),
                  "cStencilElement objects differ (pointer identity required)");
        exit(__LINE__, __FILE__, msg);
      }

      return equal;
    }
  };

  class cRhsSupportTable {
  public:
    double Coefficient;
    char *AssociatedDataPointer;

    cRhsSupportTable () {Coefficient=0.0,AssociatedDataPointer=NULL;}

    // Equality operator for class cRhsSupportTable
    bool operator==(const cRhsSupportTable& other) const {
      // Compare pointer identity for AssociatedDataPointer
      if (AssociatedDataPointer != other.AssociatedDataPointer) {
        return false;
      }

      // Compare Coefficient values for exact equality
      if (Coefficient != other.Coefficient) {
        return false;
      }

      return true;
    }

    // Inequality operator (optional but recommended)
    bool operator!=(const cRhsSupportTable& other) const {
      return !(*this == other);
    }

    // Add inside class cRhsSupportTable (public:). Requires <cstdio>.
    static bool Compare(const cRhsSupportTable& a,
                    const cRhsSupportTable& b,
                    bool break_flag = false) {
      bool equal = true;
      bool printed_header = false;

      auto hdr = [&]() {
        if (!printed_header) { std::printf("cRhsSupportTable mismatch:\n"); printed_header = true; }
      };

      // Enforce pointer identity.
      if (a.AssociatedDataPointer != b.AssociatedDataPointer) {
        equal = false; 

        if (break_flag) {std::printf("  [AssociatedDataPointer] a=%p b=%p\n",
                (const void*)a.AssociatedDataPointer, (const void*)b.AssociatedDataPointer);
	        hdr();
		exit(__LINE__,__FILE__);
	}
      }

      // Exact coefficient compare (use an epsilon if you prefer tolerance).
      if (a.Coefficient != b.Coefficient) {
        equal = false;

        if (break_flag) { std::printf("  [Coefficient] a=%.17e b=%.17e\n", a.Coefficient, b.Coefficient);
          hdr();
                exit(__LINE__,__FILE__);
        } 
      }

      if (!equal && break_flag) {
        char msg[256];
        if (break_flag) std::snprintf(msg, sizeof(msg),
                  "cRhsSupportTable objects differ (pointer identity required)");
        exit(__LINE__, __FILE__, msg);
     }

     return equal;
   }


   static bool CompareArray(const cRhsSupportTable* a,
                         const cRhsSupportTable* b,
                         int len,
                         bool break_flag = false) {
     if (len < 0) {
       std::printf("cRhsSupportTable array compare: invalid len=%d\n", len);
       if (break_flag) { char msg[64]; std::snprintf(msg,sizeof(msg),"invalid len"); exit(__LINE__, __FILE__, msg); }
    return false;
  }
  if (len == 0) return true;

  bool equal = true;
  bool printed_hdr = false;
  auto hdr = [&](){ if (!printed_hdr){ std::printf("cRhsSupportTable array mismatch:\n"); printed_hdr = true; } };

  // Track which elements of b[] have been matched.
  bool* used = new bool[len];
  memset(used, 0, sizeof(bool) * len);

  // For each a[i], find a unique equal element in b using cRhsSupportTable::Compare.
  for (int i = 0; i < len; ++i) {
    bool found = false;
    for (int j = 0; j < len; ++j) {
      if (used[j]) continue;
      if (cRhsSupportTable::Compare(a[i], b[j], false)) { // forwards break_flag
        used[j] = true;
        found = true;
        break;
      }
      // If break_flag==true and Compare found a difference, it may have exited already.
    }

    if (!found) {
      equal = false; hdr();
      std::printf("  Missing in B: ptr=%p coef=%.17e (a[%d])\n",
                  (const void*)a[i].AssociatedDataPointer, a[i].Coefficient, i);
      if (break_flag) {
        delete[] used;
        char msg[128];
        std::snprintf(msg, sizeof(msg), "Missing match in B for a[%d]", i);
        exit(__LINE__, __FILE__, msg);
      }
    }
  }

  // Report extras left unmatched in B.
  for (int j = 0; j < len; ++j) {
    if (!used[j]) {
      equal = false; hdr();
      std::printf("  Extra in B: ptr=%p coef=%.17e (b[%d])\n",
                  (const void*)b[j].AssociatedDataPointer, b[j].Coefficient, j);
    }
  }

  delete[] used;

  if (!equal && break_flag) {
    char msg[160];
    std::snprintf(msg, sizeof(msg), "cRhsSupportTable arrays differ (unordered O(N^2) compare)");
    exit(__LINE__, __FILE__, msg);
  }

  return equal;
}


  };

  class cMatrixRow {
  public:
    cMatrixRow* next;
    double Rhs;

    cStencilElement Elements[MaxStencilLength];
    cStencilElementData ElementDataTable[MaxStencilLength];
    int nNonZeroElements;

    //pointed to the beginitg of the associated data vectros used in calculating of the right-hand side vectors
    cRhsSupportTable RhsSupportTable_CornerNodes[MaxRhsSupportLength_CornerNodes];
    int RhsSupportLength_CornerNodes;

    cRhsSupportTable RhsSupportTable_CenterNodes[MaxRhsSupportLength_CenterNodes];
    int RhsSupportLength_CenterNodes;

    int i,j,k;
    int iVar;
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;
    cCornerNode* CornerNode;

    cMatrixRow() {
      next=NULL,Rhs=0.0,node=NULL,CornerNode=NULL;
      i=0,j=0,k=0,iVar=0,nNonZeroElements=0,RhsSupportLength_CornerNodes=0,RhsSupportLength_CenterNodes=0;
    }

    // -----------------------------------------------------------------
    // Equality operators
    //
    // Semantics are aligned with static Compare(a,b,...), but without
    // diagnostics or trap exits:
    //   - Pointer identity for next, node, CornerNode
    //   - Exact equality for Rhs, nNonZeroElements, i,j,k,iVar
    //   - Element-wise equality for:
    //       * Elements[0..nNonZeroElements-1]
    //       * ElementDataTable[0..nNonZeroElements-1]
    //       * RhsSupportTable_CornerNodes[0..RhsSupportLength_CornerNodes-1]
    //       * RhsSupportTable_CenterNodes[0..RhsSupportLength_CenterNodes-1]
    //
    // Requires that:
    //   - cStencilElement  implements operator==
    //   - cStencilElementData implements operator==
    //   - cRhsSupportTable implements operator==
    // -----------------------------------------------------------------
   bool operator==(const cMatrixRow& other) const {
  bool equal = true;

  // Track first mismatch info for debugging
  enum MismatchKind {
    NONE,
    POINTERS,
    SCALARS,
    STENCIL_MULTISETS,
    RHS_CORNER_MULTISETS,
    RHS_CENTER_MULTISETS
  };
  MismatchKind mismatch = NONE;
  int mismatch_index_this  = -1;
  int mismatch_index_other = -1; // currently unused, reserved if you ever want to report both sides

  // --- Pointer identity checks ---
  if (node != other.node || CornerNode != other.CornerNode) {
    equal = false;
    mismatch = POINTERS;
  }

  // --- Scalar checks ---
  if (Rhs != other.Rhs ||
      nNonZeroElements != other.nNonZeroElements ||
      i != other.i || j != other.j || k != other.k || iVar != other.iVar) {
    equal = false;
    if (mismatch == NONE) mismatch = SCALARS;
  }

  // --- Stencil payloads: (Elements, ElementDataTable) as unordered multisets ---
  if (nNonZeroElements != other.nNonZeroElements) {
    equal = false;
    if (mismatch == NONE) mismatch = STENCIL_MULTISETS;
  } else {
    bool *used = new bool[nNonZeroElements];
    for (int t = 0; t < nNonZeroElements; ++t) used[t] = false;

    for (int ia = 0; ia < nNonZeroElements; ++ia) {
      bool found  = false;
      int  foundj = -1;

      for (int jb = 0; jb < other.nNonZeroElements; ++jb) {
        if (used[jb]) continue;

        if ((Elements[ia] == other.Elements[jb]) &&
            (ElementDataTable[ia] == other.ElementDataTable[jb])) {
          used[jb] = true;
          found    = true;
          foundj   = jb;
          break;
        }
      }

      if (!found) {
        equal = false;
        if (mismatch == NONE) {
          mismatch = STENCIL_MULTISETS;
          mismatch_index_this  = ia;
          mismatch_index_other = -1;
        }
      }
    }

    delete [] used;
  }

  // --- RHS support (corner nodes) as unordered multiset ---
  if (RhsSupportLength_CornerNodes != other.RhsSupportLength_CornerNodes) {
    equal = false;
    if (mismatch == NONE) mismatch = RHS_CORNER_MULTISETS;
  } else {
    bool *used = new bool[RhsSupportLength_CornerNodes];
    for (int t = 0; t < RhsSupportLength_CornerNodes; ++t) used[t] = false;

    for (int ia = 0; ia < RhsSupportLength_CornerNodes; ++ia) {
      bool found  = false;
      int  foundj = -1;

      for (int jb = 0; jb < other.RhsSupportLength_CornerNodes; ++jb) {
        if (used[jb]) continue;
        if (RhsSupportTable_CornerNodes[ia] ==
            other.RhsSupportTable_CornerNodes[jb]) {
          used[jb] = true;
          found    = true;
          foundj   = jb;
          break;
        }
      }

      if (!found) {
        equal = false;
        if (mismatch == NONE) {
          mismatch = RHS_CORNER_MULTISETS;
          mismatch_index_this  = ia;
          mismatch_index_other = -1;
        }
      }
    }

    delete [] used;
  }

  // --- RHS support (center nodes) as unordered multiset ---
  if (RhsSupportLength_CenterNodes != other.RhsSupportLength_CenterNodes) {
    equal = false;
    if (mismatch == NONE) mismatch = RHS_CENTER_MULTISETS;
  } else {
    bool *used = new bool[RhsSupportLength_CenterNodes];
    for (int t = 0; t < RhsSupportLength_CenterNodes; ++t) used[t] = false;

    for (int ia = 0; ia < RhsSupportLength_CenterNodes; ++ia) {
      bool found  = false;
      int  foundj = -1;

      for (int jb = 0; jb < other.RhsSupportLength_CenterNodes; ++jb) {
        if (used[jb]) continue;
        if (RhsSupportTable_CenterNodes[ia] ==
            other.RhsSupportTable_CenterNodes[jb]) {
          used[jb] = true;
          found    = true;
          foundj   = jb;
          break;
        }
      }

      if (!found) {
        equal = false;
        if (mismatch == NONE) {
          mismatch = RHS_CENTER_MULTISETS;
          mismatch_index_this  = ia;
          mismatch_index_other = -1;
        }
      }
    }

    delete [] used;
  }

  // --- If not equal: dump both rows and indicate where mismatch arose ---
  if (!equal) {
    std::printf("cMatrixRow::operator== mismatch at %s:%d\n", __FILE__, __LINE__);

    // Header for this
    std::printf(
      "  this : next=%p node=%p CornerNode=%p "
      "Rhs=%.17e nNonZeroElements=%d "
      "i=%d j=%d k=%d iVar=%d "
      "RhsSupportLength_CornerNodes=%d "
      "RhsSupportLength_CenterNodes=%d\n",
      (const void*)next, (const void*)node, (const void*)CornerNode,
      Rhs, nNonZeroElements,
      i, j, k, iVar,
      RhsSupportLength_CornerNodes,
      RhsSupportLength_CenterNodes
    );

    // Header for other
    std::printf(
      "  other: next=%p node=%p CornerNode=%p "
      "Rhs=%.17e nNonZeroElements=%d "
      "i=%d j=%d k=%d iVar=%d "
      "RhsSupportLength_CornerNodes=%d "
      "RhsSupportLength_CenterNodes=%d\n",
      (const void*)other.next, (const void*)other.node, (const void*)other.CornerNode,
      other.Rhs, other.nNonZeroElements,
      other.i, other.j, other.k, other.iVar,
      other.RhsSupportLength_CornerNodes,
      other.RhsSupportLength_CenterNodes
    );

    // High-level mismatch category
    switch (mismatch) {
      case POINTERS:
        std::printf("  First mismatch category: POINTERS (next/node/CornerNode)\n");
        break;
      case SCALARS:
        std::printf("  First mismatch category: SCALARS (Rhs, nNonZeroElements or i,j,k,iVar)\n");
        break;
      case STENCIL_MULTISETS:
        std::printf("  First mismatch category: STENCIL_MULTISETS");
        if (mismatch_index_this >= 0) {
          std::printf(" – element from 'this' at idx=%d has no match in 'other'\n",
                      mismatch_index_this);
        } else {
          std::printf("\n");
        }
        break;
      case RHS_CORNER_MULTISETS:
        std::printf("  First mismatch category: RHS_CORNER_MULTISETS");
        if (mismatch_index_this >= 0) {
          std::printf(" – corner RHS support from 'this' at idx=%d has no match in 'other'\n",
                      mismatch_index_this);
        } else {
          std::printf("\n");
        }
        break;
      case RHS_CENTER_MULTISETS:
        std::printf("  First mismatch category: RHS_CENTER_MULTISETS");
        if (mismatch_index_this >= 0) {
          std::printf(" – center RHS support from 'this' at idx=%d has no match in 'other'\n",
                      mismatch_index_this);
        } else {
          std::printf("\n");
        }
        break;
      case NONE:
      default:
        break;
    }

    // --- Dump stencil payloads for both rows ---
    std::printf("  Dump of stencil (this):\n");
    for (int idx = 0; idx < nNonZeroElements; ++idx) {
      const cStencilElement&     e = Elements[idx];
      const cStencilElementData& d = ElementDataTable[idx];

      std::printf(
        "    idx=%d:\n"
        "      Element: CornerNode=%p node=%p CornerNodeID=%d "
        "ParamLen=%d SupportLen=%d\n",
        idx,
        (const void*)e.CornerNode, (const void*)e.node,
        e.CornerNodeID,
        e.MatrixElementParameterTableLength,
        e.MatrixElementSupportTableLength
      );

      std::printf("        ParamTable:");
      for (int p = 0; p < e.MatrixElementParameterTableLength; ++p) {
        std::printf(" %.17e", e.MatrixElementParameterTable[p]);
      }
      std::printf("\n");

      std::printf("        SupportTable:");
      for (int p = 0; p < e.MatrixElementSupportTableLength; ++p) {
        std::printf(" %p", e.MatrixElementSupportTable[p]);
      }
      std::printf("\n");

      std::printf(
        "      ElementData: UnknownVectorIndex=%d iVar=%d Thread=%d "
        "MatrixElementValue=%.17e rhs=%p\n",
        d.UnknownVectorIndex, d.iVar, d.Thread,
        d.MatrixElementValue, (const void*)d.rhs
      );
    }

    std::printf("  Dump of stencil (other):\n");
    for (int idx = 0; idx < other.nNonZeroElements; ++idx) {
      const cStencilElement&     e = other.Elements[idx];
      const cStencilElementData& d = other.ElementDataTable[idx];

      std::printf(
        "    idx=%d:\n"
        "      Element: CornerNode=%p node=%p CornerNodeID=%d "
        "ParamLen=%d SupportLen=%d\n",
        idx,
        (const void*)e.CornerNode, (const void*)e.node,
        e.CornerNodeID,
        e.MatrixElementParameterTableLength,
        e.MatrixElementSupportTableLength
      );

      std::printf("        ParamTable:");
      for (int p = 0; p < e.MatrixElementParameterTableLength; ++p) {
        std::printf(" %.17e", e.MatrixElementParameterTable[p]);
      }
      std::printf("\n");

      std::printf("        SupportTable:");
      for (int p = 0; p < e.MatrixElementSupportTableLength; ++p) {
        std::printf(" %p", e.MatrixElementSupportTable[p]);
      }
      std::printf("\n");

      std::printf(
        "      ElementData: UnknownVectorIndex=%d iVar=%d Thread=%d "
        "MatrixElementValue=%.17e rhs=%p\n",
        d.UnknownVectorIndex, d.iVar, d.Thread,
        d.MatrixElementValue, (const void*)d.rhs
      );
    }

    // --- Dump RHS support tables for both rows ---
    std::printf("  RHS support (corner, this):\n");
    for (int n = 0; n < RhsSupportLength_CornerNodes; ++n) {
      std::printf("    idx=%d: Coefficient=%.17e ptr=%p\n",
                  n,
                  RhsSupportTable_CornerNodes[n].Coefficient,
                  (const void*)RhsSupportTable_CornerNodes[n].AssociatedDataPointer);
    }
    std::printf("  RHS support (corner, other):\n");
    for (int n = 0; n < other.RhsSupportLength_CornerNodes; ++n) {
      std::printf("    idx=%d: Coefficient=%.17e ptr=%p\n",
                  n,
                  other.RhsSupportTable_CornerNodes[n].Coefficient,
                  (const void*)other.RhsSupportTable_CornerNodes[n].AssociatedDataPointer);
    }

    std::printf("  RHS support (center, this):\n");
    for (int n = 0; n < RhsSupportLength_CenterNodes; ++n) {
      std::printf("    idx=%d: Coefficient=%.17e ptr=%p\n",
                  n,
                  RhsSupportTable_CenterNodes[n].Coefficient,
                  (const void*)RhsSupportTable_CenterNodes[n].AssociatedDataPointer);
    }
    std::printf("  RHS support (center, other):\n");
    for (int n = 0; n < other.RhsSupportLength_CenterNodes; ++n) {
      std::printf("    idx=%d: Coefficient=%.17e ptr=%p\n",
                  n,
                  other.RhsSupportTable_CenterNodes[n].Coefficient,
                  (const void*)other.RhsSupportTable_CenterNodes[n].AssociatedDataPointer);
    }
  }

  return equal;
}
 
    
    
    bool operator!=(const cMatrixRow& other) const {
      return !(*this == other);
    }
  };

  // Add inside class cMatrixRow (public:). Requires <cstdio>.
  // Assumes cStencilElement::Compare(...), cStencilElementData::Compare(...),
  // and cRhsSupportTable::Compare(...) are defined as in earlier steps.
  static bool Compare(const cMatrixRow& a,
                    const cMatrixRow& b,
                    bool break_flag = false) {
    bool equal = true;
    bool printed_header = false;
    auto hdr = [&](){ if (!printed_header) { std::printf("cMatrixRow mismatch:\n"); printed_header = true; } };

    // --- Pointer identity checks ---
    if (a.next != b.next) { equal = false; hdr();
      std::printf("  [next ptr] a=%p b=%p\n", (const void*)a.next, (const void*)b.next);
    }
    if (a.node != b.node) { equal = false; hdr();
      std::printf("  [node ptr] a=%p b=%p\n", (const void*)a.node, (const void*)b.node);
    }
    if (a.CornerNode != b.CornerNode) { equal = false; hdr();
      std::printf("  [CornerNode ptr] a=%p b=%p\n", (const void*)a.CornerNode, (const void*)b.CornerNode);
    }

    // --- Scalars ---
    if (a.Rhs != b.Rhs) { equal = false; hdr();
      std::printf("  [Rhs] a=%.17e b=%.17e\n", a.Rhs, b.Rhs);
    }
    if (a.nNonZeroElements != b.nNonZeroElements) { equal = false; hdr();
      std::printf("  [nNonZeroElements] a=%d b=%d\n", a.nNonZeroElements, b.nNonZeroElements);
    }
    if (a.i != b.i) { equal = false; hdr(); std::printf("  [i] a=%d b=%d\n", a.i, b.i); }
    if (a.j != b.j) { equal = false; hdr(); std::printf("  [j] a=%d b=%d\n", a.j, b.j); }
    if (a.k != b.k) { equal = false; hdr(); std::printf("  [k] a=%d b=%d\n", a.k, b.k); }
    if (a.iVar != b.iVar) { equal = false; hdr(); std::printf("  [iVar] a=%d b=%d\n", a.iVar, b.iVar); }

    // --- Stencil payloads: compare up to min(nNonZeroElements) ---
    {
      const int nnz_min = (a.nNonZeroElements < b.nNonZeroElements) ? a.nNonZeroElements : b.nNonZeroElements;
      for (int idx = 0; idx < nnz_min; ++idx) {
        // Elements (enforce pointer identity inside the element comparer)
        if (!cStencilElement::Compare(a.Elements[idx], b.Elements[idx], /*break_flag*/ false)) {
          equal = false; hdr();
          std::printf("  [Elements[%d]] differ\n", idx);
        }
        // ElementDataTable (compares indices, value, and rhs pointer identity)
        if (!cStencilElementData::Compare(a.ElementDataTable[idx], b.ElementDataTable[idx], /*break_flag*/ false)) {
          equal = false; hdr();
          std::printf("  [ElementDataTable[%d]] differ\n", idx);
        }
      }
      // If counts differ, flag the tail
      for (int idx = nnz_min; idx < a.nNonZeroElements; ++idx) { equal = false; hdr();
        std::printf("  [Elements[%d]/ElementDataTable[%d]] present only in 'a'\n", idx, idx);
      }
      for (int idx = nnz_min; idx < b.nNonZeroElements; ++idx) { equal = false; hdr();
        std::printf("  [Elements[%d]/ElementDataTable[%d]] present only in 'b'\n", idx, idx);
      }
    }

    // --- RHS support (corner nodes) ---
    if (a.RhsSupportLength_CornerNodes != b.RhsSupportLength_CornerNodes) { equal = false; hdr();
      std::printf("  [RhsSupportLength_CornerNodes] a=%d b=%d\n",
                a.RhsSupportLength_CornerNodes, b.RhsSupportLength_CornerNodes);
    }
    {
      const int len_min = (a.RhsSupportLength_CornerNodes < b.RhsSupportLength_CornerNodes)
                        ? a.RhsSupportLength_CornerNodes : b.RhsSupportLength_CornerNodes;
      for (int i = 0; i < len_min; ++i) {
        if (!cRhsSupportTable::Compare(a.RhsSupportTable_CornerNodes[i],
                                     b.RhsSupportTable_CornerNodes[i],
                                     /*break_flag*/ false)) {
          equal = false; hdr();
          std::printf("  [RhsSupportTable_CornerNodes[%d]] differ\n", i);
        }
      }
      for (int i = len_min; i < a.RhsSupportLength_CornerNodes; ++i) { equal = false; hdr();
        std::printf("  [RhsSupportTable_CornerNodes[%d]] present only in 'a'\n", i);
      }
      for (int i = len_min; i < b.RhsSupportLength_CornerNodes; ++i) { equal = false; hdr();
        std::printf("  [RhsSupportTable_CornerNodes[%d]] present only in 'b'\n", i);
      }
    }

    // --- RHS support (center nodes) ---
    if (a.RhsSupportLength_CenterNodes != b.RhsSupportLength_CenterNodes) { equal = false; hdr();
      std::printf("  [RhsSupportLength_CenterNodes] a=%d b=%d\n",
                a.RhsSupportLength_CenterNodes, b.RhsSupportLength_CenterNodes);
    }
    {
      const int len_min = (a.RhsSupportLength_CenterNodes < b.RhsSupportLength_CenterNodes)
                        ? a.RhsSupportLength_CenterNodes : b.RhsSupportLength_CenterNodes;
      for (int i = 0; i < len_min; ++i) {
        if (!cRhsSupportTable::Compare(a.RhsSupportTable_CenterNodes[i],
                                     b.RhsSupportTable_CenterNodes[i],
                                     /*break_flag*/ false)) {
          equal = false; hdr();
          std::printf("  [RhsSupportTable_CenterNodes[%d]] differ\n", i);
        }
      }
      for (int i = len_min; i < a.RhsSupportLength_CenterNodes; ++i) { equal = false; hdr();
        std::printf("  [RhsSupportTable_CenterNodes[%d]] present only in 'a'\n", i);
      }
      for (int i = len_min; i < b.RhsSupportLength_CenterNodes; ++i) { equal = false; hdr();
        std::printf("  [RhsSupportTable_CenterNodes[%d]] present only in 'b'\n", i);
      }
    }

    if (!equal && break_flag) {
      char msg[256];
      std::snprintf(msg, sizeof(msg), "cMatrixRow objects differ (see printf log above)");
      exit(__LINE__, __FILE__, msg);
    }

    return equal;
  }


  static bool CompareArray(const cMatrixRow* A,
                         const cMatrixRow* B,
                         int len,
                         bool break_flag = false) {
    if (len < 0) {
      std::printf("cMatrixRow array compare: invalid len=%d\n", len);
      if (break_flag) { char msg[64]; std::snprintf(msg,sizeof(msg),"invalid len"); exit(__LINE__, __FILE__, msg); }
      return false;
    }
    if (len == 0) return true;

    bool equal = true;
    bool printed_hdr = false;
    auto hdr = [&](){ if (!printed_hdr){ std::printf("cMatrixRow array mismatch:\n"); printed_hdr = true; } };

    // Track which B[j] have been matched already
    bool* used = new bool[len];
    memset(used, 0, sizeof(bool)*len);

    // For each A[i], find a unique equal row in B using cMatrixRow::Compare
    for (int i = 0; i < len; ++i) {
      bool found = false;
      for (int j = 0; j < len; ++j) {
        if (used[j]) continue;
        if (cMatrixRow::Compare(A[i], B[j], false)) {
          used[j] = true;
          found = true;
          break;
        }
        // If break_flag==true and Compare found a difference, it may have exited already.
      }

      if (!found) {
        equal = false; hdr();
        std::printf("  Missing in B: key=(%d,%d,%d|%d)  (A[%d])\n",
                  A[i].i, A[i].j, A[i].k, A[i].iVar, i);
        if (break_flag) {
          delete[] used;
          char msg[128];
          std::snprintf(msg, sizeof(msg), "Missing match in B for A[%d] (key=(%d,%d,%d|%d))",
                      i, A[i].i, A[i].j, A[i].k, A[i].iVar);
          exit(__LINE__, __FILE__, msg);
        }
      }
    }

    // Any unused entries in B are extras
    for (int j = 0; j < len; ++j) {
      if (!used[j]) {
        equal = false; hdr();
        std::printf("  Extra in B: key=(%d,%d,%d|%d)  (B[%d])\n",
                  B[j].i, B[j].j, B[j].k, B[j].iVar, j);
      }
    }

    delete[] used;

    if (!equal && break_flag) {
      char msg[160];
      std::snprintf(msg, sizeof(msg), "cMatrixRow arrays differ (unordered O(N^2) compare)");
      exit(__LINE__, __FILE__, msg);
    }

    return equal;
  }


  cStackManaged<cMatrixRow> MatrixRowStack;

  //Table of the rows local to the current MPI process
  cMatrixRow *MatrixRowListFirst,*MatrixRowListLast;
  cMatrixRow **MatrixRowTable;
  int MatrixRowTableLength;

  //exchange buffers
  int* RecvExchangeBufferLength;  //the number of elementf in the recv list
  int* SendExchangeBufferLength;
  int **SendExchangeBufferElementIndex;

  //add new row to the matrix
  //for that a user function is called that returns stenal (elements of the matrix), corresponding rhs, and the number of the elements in the amtrix' line
  class cMatrixRowNonZeroElementTable {
  public:
    int i,j,k; //coordintes of the corner block used in the equation
    int iVar; //the index of the variable used in the stencil
    double MatrixElementValue;
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *Node;

    //the following are needed only in case when the periodic bounday conditinos are imposed
    bool BoundaryNodeFlag;
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *OriginalNode;

    //the user-model parameters neede for recalculating of the matrix coefficients
    double MatrixElementParameterTable[MaxMatrixElementParameterTableLength];
    int MatrixElementParameterTableLength;

    //the user-defined support used for recalculating the matrix coefficients
    void *MatrixElementSupportTable[MaxMatrixElementSupportTableLength];
    int MatrixElementSupportTableLength;

    cMatrixRowNonZeroElementTable() {
      i=0,j=0,k=0,iVar=0,MatrixElementValue=0.0,Node=NULL;
      BoundaryNodeFlag=false,OriginalNode=NULL;

      MatrixElementParameterTableLength=0,MatrixElementSupportTableLength=0;

      for (int i=0;i<MaxMatrixElementParameterTableLength;i++) MatrixElementParameterTable[i]=0.0;
      for (int i=0;i<MaxMatrixElementSupportTableLength;i++) MatrixElementSupportTable[i]=NULL;
    }

    // -----------------------------------------------------------------
    // Equality operators
    //
    // Two cMatrixRowNonZeroElementTable entries are considered equal iff:
    //   - i,j,k,iVar match exactly
    //   - MatrixElementValue matches exactly (no tolerance)
    //   - Node, BoundaryNodeFlag, OriginalNode match
    //   - MatrixElementParameterTableLength matches and all entries
    //     [0 .. length-1] are exactly equal
    //   - MatrixElementSupportTableLength matches and all entries
    //     [0 .. length-1] are pointer-identical
    // -----------------------------------------------------------------
    bool operator==(const cMatrixRowNonZeroElementTable& other) const {
      // Indices and value
      if (i   != other.i)   return false;
      if (j   != other.j)   return false;
      if (k   != other.k)   return false;
      if (iVar != other.iVar) return false;
      if (MatrixElementValue != other.MatrixElementValue) return false;

      // Node and periodic info
      if (Node           != other.Node)           return false;
      if (BoundaryNodeFlag != other.BoundaryNodeFlag) return false;
      if (OriginalNode   != other.OriginalNode)   return false;

      // Parameter table
      if (MatrixElementParameterTableLength !=
          other.MatrixElementParameterTableLength) {
        return false;
      }

      for (int n = 0; n < MatrixElementParameterTableLength; ++n) {
        if (MatrixElementParameterTable[n] !=
          other.MatrixElementParameterTable[n]) {
          return false;
        }
      }

      // Support table (pointer identity)
      if (MatrixElementSupportTableLength !=
        other.MatrixElementSupportTableLength) {
        return false;
      }
      
      for (int n = 0; n < MatrixElementSupportTableLength; ++n) {
        if (MatrixElementSupportTable[n] !=
          other.MatrixElementSupportTable[n]) {
          return false;
        }
      }

      return true;
    }

    bool operator!=(const cMatrixRowNonZeroElementTable& other) const {
      return !(*this == other);
    }
  };

  cMatrixRowNonZeroElementTable MatrixRowNonZeroElementTable[MaxStencilLength];

  //provcess the right boundary
  int nVirtualBlocks,nRealBlocks;
  void ProcessRightDomainBoundary(int *RecvDataPointCounter,void(*f)(int i,int j,int k,int iVar,cMatrixRowNonZeroElementTable* Set,int& NonZeroElementsFound,double& Rhs,cRhsSupportTable* RhsSupportTable_CornerNodes,int &RhsSupportLength_CornerNodes,cRhsSupportTable* RhsSupportTable_CenterNodes,int &RhsSupportLength_CenterNodes,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node));


  //void reset the data of the obsect to the default state
  _TARGET_HOST_ _TARGET_DEVICE_
  void DeleteDataBuffers();

  _TARGET_HOST_ _TARGET_DEVICE_
  void Reset(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode);

  _TARGET_HOST_ _TARGET_DEVICE_
  void Reset();

  //reset indexing of the nodes
  _TARGET_HOST_ _TARGET_DEVICE_
  void ResetUnknownVectorIndex(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node);

  //build the matrix
  void(*GetStencil)(int,int,int,int,cMatrixRowNonZeroElementTable*,int&,double&,cRhsSupportTable*,int&,cRhsSupportTable*,int&,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*); 

  _TARGET_HOST_ 
  void BuildMatrix(void(*f)(int i,int j,int k,int iVar,cMatrixRowNonZeroElementTable* Set,int& NonZeroElementsFound,double& Rhs,cRhsSupportTable* RhsSupportTable_CornerNodes,int &RhsSupportLength_CornerNodes,cRhsSupportTable* RhsSupportTable_CenterNodes,int &RhsSupportLength_CenterNodes,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node));

  //re-build the matrix after the domain re-decomposition operation
  void RebuildMatrix() {
    if (GetStencil!=NULL) BuildMatrix(GetStencil);
  }

  //constructor
  cLinearSystemCornerNode() {
    MatrixRowListFirst=NULL,MatrixRowListLast=NULL;
    MatrixRowTable=NULL,MatrixRowTableLength=0;

    //counter of the matrix modifications
    nMartixModifications=0;

    //exchange buffers
    RecvExchangeBufferLength=NULL;
    SendExchangeBufferLength=NULL,SendExchangeBufferElementIndex=NULL;

    //partial data of the linear system to be solved
    SubdomainPartialRHS=NULL,SubdomainPartialUnknownsVector=NULL;
    SubdomainPartialUnknownsVectorLength=0;

    nVirtualBlocks=0,nRealBlocks=0;
    GetStencil=NULL;
  }

  //exchange the data
  void ExchangeIntermediateUnknownsData(double *x,double** &);

  //matrix/vector multiplication
  _TARGET_HOST_
  void MultiplyVector(double *p,double *x,int length);

  //call the linear system solver, and unpack the solution afterward
  void Solve(void (*fInitialUnknownValues)(double* x,cCornerNode* CornerNode),void (*fUnpackSolution)(double* x,cCornerNode* CornerNode),double Tol,int nMaxIter,
      int (*fPackBlockData)(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>** NodeTable,int NodeTableLength,unsigned long int* NodeDataLength,unsigned char* BlockCenterNodeMask,unsigned char* BlockCornerNodeMask,char* SendDataBuffer),
      int (*fUnpackBlockData)(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>** NodeTable,int NodeTableLength,unsigned char* BlockCenterNodeMask,unsigned char* BlockCornerNodeMask,char* RecvDataBuffer));

  //update the RHS vector
  void UpdateRhs(double (*fSetRhs)(int,cRhsSupportTable*,int,cRhsSupportTable*,int));

  //update non-zero coefficients of the matrix
  void UpdateMatrixNonZeroCoefficients(void (*UpdateMatrixRow)(cMatrixRow*));

  //destructor
  ~cLinearSystemCornerNode() {
    if (MatrixRowTable!=NULL) delete [] MatrixRowTable;
  }

  //calculate signature of the matrix
  void GetSignature(int nline,const char* fname) {
    CRC32 Signature,SignatureRhs;
    cMatrixRow* row;
    int cnt,iElementMax,iElement;
    cStencilElementData *data,*ElementDataTable;

    for (row=MatrixRowListFirst,cnt=0;row!=NULL;row=row->next,cnt++) {
      Signature.add(cnt);

      SignatureRhs.add(cnt);
      SignatureRhs.add(row->Rhs);

      iElementMax=row->nNonZeroElements;
      ElementDataTable=row->ElementDataTable;

      for (iElement=0;iElement<iElementMax;iElement++) {
        data=ElementDataTable+iElement;

        Signature.add(data->MatrixElementValue);
        Signature.add(data->Thread);
        Signature.add(data->iVar);
        Signature.add(data->UnknownVectorIndex);
      }
    }

    Signature.PrintChecksum(nline,fname);
    SignatureRhs.PrintChecksum(nline,fname);
  }

  // Strict structural equality for the linear system.
  //
  // Notes:
  //  - Uses cMatrixRow::operator== for each row.
  //  - Compares SubdomainPartialRHS element-wise when both pointers are non-null.
  //  - Pointer-valued members that encode ownership/layout are compared by identity,
  //    consistent with the nested comparison helpers in this file.
  bool operator==(const cLinearSystemCornerNode& other) const {
    // Basic configuration
    if (SubdomainPartialUnknownsVectorLength != other.SubdomainPartialUnknownsVectorLength)
      return false;

    if (nMartixModifications != other.nMartixModifications)
      return false;

    // Partial RHS / solution buffers: compare identity first
    if (SubdomainPartialRHS != other.SubdomainPartialRHS)
      return false;

    if (SubdomainPartialUnknownsVector != other.SubdomainPartialUnknownsVector)
      return false;

    // Local row table configuration
    if (MatrixRowTableLength != other.MatrixRowTableLength)
      return false;

    if (MatrixRowListFirst != other.MatrixRowListFirst)
      return false;

    if (MatrixRowListLast != other.MatrixRowListLast)
      return false;

    // Exchange buffers – pointer identity
    if (RecvExchangeBufferLength != other.RecvExchangeBufferLength)
      return false;

    if (SendExchangeBufferLength != other.SendExchangeBufferLength)
      return false;

    if (SendExchangeBufferElementIndex != other.SendExchangeBufferElementIndex)
      return false;

    // Matrix rows: use cMatrixRow::operator==
    for (int i = 0; i < MatrixRowTableLength; ++i) {
      cMatrixRow* rowA = MatrixRowTable[i];
      cMatrixRow* rowB = other.MatrixRowTable[i];

      if (rowA == nullptr || rowB == nullptr) {
        if (rowA != rowB) return false; // one null, one non-null
        continue;                       // both null -> ok
      }

      if (!(*rowA == *rowB)) {
        return false;
      }
    }

    // RHS vector values (if allocated in both objects and length is consistent)
    if (SubdomainPartialRHS && other.SubdomainPartialRHS &&
        SubdomainPartialUnknownsVectorLength == other.SubdomainPartialUnknownsVectorLength) {
      for (int i = 0; i < SubdomainPartialUnknownsVectorLength; ++i) {
        if (SubdomainPartialRHS[i] != other.SubdomainPartialRHS[i]) {
          return false;
        }
      }
    }

    return true;
  }

  bool operator!=(const cLinearSystemCornerNode& other) const {
    return !(*this == other);
  }

  // -----------------------------------------------------------------
  // Debug-style comparison with optional hard exit.
  //
  // If exit_flag == true and any difference is detected, this prints
  // a short log to stdout and then calls:
  //     exit(__LINE__, __FILE__, msg);
  //
  // Return value:
  //   - true  : objects are equal (under the same semantics as operator==)
  //   - false : objects differ
  // -----------------------------------------------------------------
// -----------------------------------------------------------------
// Debug-style comparison for entire linear systems (pointer version).
//
// Parameters:
//   a, b      - pointers to the systems to compare
//   exit_flag - if true, call exit(__LINE__, __FILE__, msg) on mismatch
//
// Return:
//   true  - systems are equal under the same semantics as operator==
//   false - mismatch found (unless exit_flag==true, in which case
//           the function will not return and will call exit()).
// -----------------------------------------------------------------
static bool Compare(cLinearSystemCornerNode* a,
                    cLinearSystemCornerNode* b,
                    bool exit_flag = false) {
  bool equal = true;
  bool printed_header = false;

  auto hdr = [&]() {
    if (!printed_header) {
      std::printf("cLinearSystemCornerNode mismatch:\n");
      printed_header = true;
    }
  };

  // ---------- Null pointer handling ----------
  if (a == nullptr || b == nullptr) {
    if (a == b) {
      // both nullptr -> treat as equal
      return true;
    }

    equal = false;
    hdr();
    std::printf("  [pointer] one system is nullptr, the other is not (a=%p b=%p)\n",
                (const void*)a, (const void*)b);

    if (exit_flag) {
      char msg[256];
      std::snprintf(msg, sizeof(msg),
                    "cLinearSystemCornerNode::Compare - nullptr vs non-null");
      exit(__LINE__, __FILE__, msg);
    }

    return false;
  }

  // For brevity, use references internally
  const cLinearSystemCornerNode& A = *a;
  const cLinearSystemCornerNode& B = *b;

  // ---------- Basic configuration ----------
  if (A.SubdomainPartialUnknownsVectorLength !=
      B.SubdomainPartialUnknownsVectorLength) {
    equal = false; hdr();
    std::printf("  [SubdomainPartialUnknownsVectorLength] a=%d b=%d\n",
                A.SubdomainPartialUnknownsVectorLength,
                B.SubdomainPartialUnknownsVectorLength);
  }

  if (A.nMartixModifications != B.nMartixModifications) {
    equal = false; hdr();
    std::printf("  [nMartixModifications] a=%d b=%d\n",
                A.nMartixModifications, B.nMartixModifications);
  }

  // ---------- Matrix rows (structural compare) ----------
  {
    const int len_min =
        (A.MatrixRowTableLength < B.MatrixRowTableLength)
            ? A.MatrixRowTableLength
            : B.MatrixRowTableLength;

    for (int i = 0; i < len_min; ++i) {
      typename cLinearSystemCornerNode::cMatrixRow* rowA = A.MatrixRowTable[i];
      typename cLinearSystemCornerNode::cMatrixRow* rowB = B.MatrixRowTable[i];

      if (rowA == nullptr || rowB == nullptr) {
        if (rowA != rowB) {
          equal = false; hdr();
          std::printf("  [MatrixRowTable[%d]] one null, one non-null (a=%p b=%p)\n",
                      i, (const void*)rowA, (const void*)rowB);
        }
        continue;
      }

      // Use the row-level comparator for detailed diagnostics
      if (*rowA!=*rowB) {
        equal = false; hdr();
        std::printf("  [MatrixRowTable[%d]] rows differ (see cMatrixRow mismatch above)\n", i);
      }
    }

    for (int i = len_min; i < A.MatrixRowTableLength; ++i) {
      equal = false; hdr();
      std::printf("  [MatrixRowTable[%d]] present only in 'a'\n", i);
    }
    for (int i = len_min; i < B.MatrixRowTableLength; ++i) {
      equal = false; hdr();
      std::printf("  [MatrixRowTable[%d]] present only in 'b'\n", i);
    }
  }

  // ---------- RHS vector values (element-wise compare) ----------
  if (A.SubdomainPartialRHS && B.SubdomainPartialRHS &&
      A.SubdomainPartialUnknownsVectorLength ==
      B.SubdomainPartialUnknownsVectorLength) {

    for (int i = 0; i < A.SubdomainPartialUnknownsVectorLength; ++i) {
      if (A.SubdomainPartialRHS[i] != B.SubdomainPartialRHS[i]) {
        equal = false; hdr();
        std::printf("  [SubdomainPartialRHS[%d]] a=%.17e b=%.17e\n",
                    i, A.SubdomainPartialRHS[i], B.SubdomainPartialRHS[i]);
        break; // one mismatch is enough for diagnostics
      }
    }
  }

  if (!equal && exit_flag) {
    char msg[256];
    std::snprintf(msg, sizeof(msg),
                  "cLinearSystemCornerNode objects differ (see printf log above)");
    exit(__LINE__, __FILE__, msg);
  }

  return equal;
}


};

template <class cCornerNode, int NodeUnknownVariableVectorLength,int MaxStencilLength,
int MaxRhsSupportLength_CornerNodes,int MaxRhsSupportLength_CenterNodes,
int MaxMatrixElementParameterTableLength,int MaxMatrixElementSupportTableLength>
_TARGET_HOST_ _TARGET_DEVICE_
void cLinearSystemCornerNode<cCornerNode, NodeUnknownVariableVectorLength,MaxStencilLength,MaxRhsSupportLength_CornerNodes,MaxRhsSupportLength_CenterNodes,MaxMatrixElementParameterTableLength,MaxMatrixElementSupportTableLength>::ResetUnknownVectorIndex(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  if (node->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
    PIC::Mesh::cDataBlockAMR *block=NULL;
    PIC::Mesh::cDataCornerNode *CornerNode=NULL;

    if ((block=node->block)!=NULL) {
      for (int i=0;i<_BLOCK_CELLS_X_;i++) for (int j=0;j<_BLOCK_CELLS_Y_;j++) for (int k=0;k<_BLOCK_CELLS_Z_;k++) {
        if ((CornerNode=block->GetCornerNode(PIC::Mesh::mesh->getCornerNodeLocalNumber(i,j,k)))!=NULL) {
          CornerNode->LinearSolverUnknownVectorIndex=-1;
        }
      }
    }

  }
  else for (int i=0;i<(1<<DIM);i++) if (node->downNode[i]!=NULL) ResetUnknownVectorIndex(node->downNode[i]);
}

template <class cCornerNode, int NodeUnknownVariableVectorLength,int MaxStencilLength,
int MaxRhsSupportLength_CornerNodes,int MaxRhsSupportLength_CenterNodes,
int MaxMatrixElementParameterTableLength,int MaxMatrixElementSupportTableLength>
void cLinearSystemCornerNode<cCornerNode, NodeUnknownVariableVectorLength,MaxStencilLength,MaxRhsSupportLength_CornerNodes,MaxRhsSupportLength_CenterNodes,MaxMatrixElementParameterTableLength,MaxMatrixElementSupportTableLength>::ProcessRightDomainBoundary(int* RecvDataPointCounter,void(*f)(int i,int j,int k,int iVar,cMatrixRowNonZeroElementTable* Set,int& NonZeroElementsFound,double& Rhs,cRhsSupportTable* RhsSupportTable_CornerNodes,int &RhsSupportLength_CornerNodes,cRhsSupportTable* RhsSupportTable_CenterNodes,int &RhsSupportLength_CenterNodes,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node)) {
  int i,j,k,iface,iblock;
  cMatrixRow* NewRow;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;
  cStencilElement* el;
  cStencilElementData* el_data;
//  int Index=PIC::DomainBlockDecomposition::nLocalBlocks*_BLOCK_CELLS_Z_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_X_;
  int kMax=_BLOCK_CELLS_Z_,jMax=_BLOCK_CELLS_Y_,iMax=_BLOCK_CELLS_X_;

//  int Debug=0;

  struct cOffset {
    int di,dj,dk;
  };


  nVirtualBlocks=0;

//  return ;

  cOffset xOffset[1]={{_BLOCK_CELLS_X_,0,0}};
  cOffset yOffset[1]={{0,_BLOCK_CELLS_Y_,0}};
  cOffset zOffset[1]={{0,0,_BLOCK_CELLS_Z_}};

  cOffset xyOffset[3]={{_BLOCK_CELLS_X_,0,0},{0,_BLOCK_CELLS_Y_,0},{_BLOCK_CELLS_X_,_BLOCK_CELLS_Y_,0}};
  cOffset xzOffset[3]={{_BLOCK_CELLS_X_,0,0},{0,0,_BLOCK_CELLS_Z_},{_BLOCK_CELLS_X_,0,_BLOCK_CELLS_Z_}};
  cOffset yzOffset[3]={{0,_BLOCK_CELLS_Y_,0},{0,0,_BLOCK_CELLS_Z_},{0,_BLOCK_CELLS_Y_,_BLOCK_CELLS_Z_}};

  cOffset xyzOffset[7]={{_BLOCK_CELLS_X_,0,0},{0,_BLOCK_CELLS_Y_,0},{0,0,_BLOCK_CELLS_Z_},
      {_BLOCK_CELLS_X_,_BLOCK_CELLS_Y_,0},{_BLOCK_CELLS_X_,0,_BLOCK_CELLS_Z_},{0,_BLOCK_CELLS_Y_,_BLOCK_CELLS_Z_},
      {_BLOCK_CELLS_X_,_BLOCK_CELLS_Y_,_BLOCK_CELLS_Z_}};

  for (int nLocalNode=0;nLocalNode<PIC::DomainBlockDecomposition::nLocalBlocks;nLocalNode++) {
    bool flag=false;

    node=PIC::DomainBlockDecomposition::BlockTable[nLocalNode];

    if (node->Thread==PIC::ThisThread) for (iface=1;iface<6;iface+=2) if (node->GetNeibFace(iface,0,0,PIC::Mesh::mesh)==NULL) {
      flag=true;
      break;
    }

    if (flag==false) continue;

/*Debug++;

//debug -> count the number of rows and blocks
int ntotbl=0;
for ( cMatrixRow* Row=MatrixRowListFirst;Row!=NULL;Row=Row->next) ntotbl++;*/

    //the block has at least one face at the 'right' boundary of the domain
    //determine the combination of the faces that are at the boundary
    int xBoundaryFace=0,yBoundaryFace=0,zBoundaryFace=0;

    if (node->GetNeibFace(1,0,0,PIC::Mesh::mesh)==NULL) xBoundaryFace=10;
    if (node->GetNeibFace(3,0,0,PIC::Mesh::mesh)==NULL) yBoundaryFace=100;
    if (node->GetNeibFace(5,0,0,PIC::Mesh::mesh)==NULL) zBoundaryFace=1000;

    //build "virtual" blocks for the combination of the boundaries
    cOffset *OffsetTable;
    int OffsetTableLength;

    switch (xBoundaryFace+yBoundaryFace+zBoundaryFace) {
     case 10+0+0:
      OffsetTable=xOffset;
      OffsetTableLength=1;
      break;

     case 0+100+0:
      OffsetTable=yOffset;
      OffsetTableLength=1;
      break;

     case 0+0+1000:
      OffsetTable=zOffset;
      OffsetTableLength=1;
      break;

     case 10+100+0:
      OffsetTable=xyOffset;
      OffsetTableLength=3;
      break;

     case 10+0+1000:
      OffsetTable=xzOffset;
      OffsetTableLength=3;
      break;

     case 0+100+1000:
      OffsetTable=yzOffset;
      OffsetTableLength=3;
      break;

     case 10+100+1000:
      OffsetTable=xyzOffset;
      OffsetTableLength=7;
      break;
    }

    nVirtualBlocks+=OffsetTableLength;

    //loop through all 'virtual' blocks
    for (iblock=0;iblock<OffsetTableLength;iblock++) {
      //loop through all 'internal' corner nodes of the block
      for (k=0;k<_BLOCK_CELLS_Z_;k++) for (j=0;j<_BLOCK_CELLS_Y_;j++) for (i=0;i<_BLOCK_CELLS_X_;i++) for (int iVar=0;iVar<NodeUnknownVariableVectorLength;iVar++) {
        if ( (i+OffsetTable[iblock].di<=_BLOCK_CELLS_X_) && (j+OffsetTable[iblock].dj<=_BLOCK_CELLS_Y_) && (k+OffsetTable[iblock].dk<=_BLOCK_CELLS_Z_) ) {
          //the point is located at the boundary of the domain and the user-defined funtion is needed to be called to get the interpolation stencil
          int NonZeroElementsFound;
          double rhs;

          NewRow=MatrixRowStack.newElement();

          //call the user defined function to determine the non-zero elements of the matrix
          NewRow->RhsSupportLength_CornerNodes=0;
          NewRow->RhsSupportLength_CenterNodes=0;
          rhs=0.0;

          for (int ii=0;ii<MaxStencilLength;ii++) MatrixRowNonZeroElementTable[ii].MatrixElementParameterTableLength=0,MatrixRowNonZeroElementTable[ii].MatrixElementSupportTableLength=0;

          f(i+OffsetTable[iblock].di,j+OffsetTable[iblock].dj,k+OffsetTable[iblock].dk,iVar,MatrixRowNonZeroElementTable,NonZeroElementsFound,rhs,NewRow->RhsSupportTable_CornerNodes,NewRow->RhsSupportLength_CornerNodes,NewRow->RhsSupportTable_CenterNodes,NewRow->RhsSupportLength_CenterNodes,node);
          if (NonZeroElementsFound>MaxStencilLength) exit(__LINE__,__FILE__,"Error: NonZeroElementsFound>=nMaxMatrixNonzeroElement; Need to increase the value of nMaxMatrixNonzeroElement");

          //populate the new row
          NewRow->i=i+OffsetTable[iblock].di,NewRow->j=j+OffsetTable[iblock].dj,NewRow->k=k+OffsetTable[iblock].dk;
          NewRow->iVar=iVar;
          NewRow->node=node;
          NewRow->Rhs=rhs;
          NewRow->CornerNode=node->block->GetCornerNode(PIC::Mesh::mesh->getCornerNodeLocalNumber(NewRow->i,NewRow->j,NewRow->k));
          NewRow->nNonZeroElements=NonZeroElementsFound;

          //add the new row to the matrix
          if (MatrixRowListFirst==NULL) {
            MatrixRowListLast=NewRow;
            MatrixRowListFirst=NewRow;
            NewRow->next=NULL;
          }
          else {
            MatrixRowListLast->next=NewRow;
            MatrixRowListLast=NewRow;
            NewRow->next=NULL;
          }

          //add to the row non-zero elements
          for (int ii=0;ii<NonZeroElementsFound;ii++) {
            MatrixRowNonZeroElementTable[ii].Node=node;

            if ((MatrixRowNonZeroElementTable[ii].i>=iMax) && (MatrixRowNonZeroElementTable[ii].Node->GetNeibFace(1,0,0,PIC::Mesh::mesh)!=NULL)) {
              MatrixRowNonZeroElementTable[ii].i-=_BLOCK_CELLS_X_;
              MatrixRowNonZeroElementTable[ii].Node=MatrixRowNonZeroElementTable[ii].Node->GetNeibFace(1,0,0,PIC::Mesh::mesh);
            }
            else if (MatrixRowNonZeroElementTable[ii].i<0) {
              MatrixRowNonZeroElementTable[ii].i+=_BLOCK_CELLS_X_;
              MatrixRowNonZeroElementTable[ii].Node=MatrixRowNonZeroElementTable[ii].Node->GetNeibFace(0,0,0,PIC::Mesh::mesh);
            }

            if ((MatrixRowNonZeroElementTable[ii].j>=jMax) && (MatrixRowNonZeroElementTable[ii].Node->GetNeibFace(3,0,0,PIC::Mesh::mesh)!=NULL)) {
              MatrixRowNonZeroElementTable[ii].j-=_BLOCK_CELLS_Y_;
              MatrixRowNonZeroElementTable[ii].Node=MatrixRowNonZeroElementTable[ii].Node->GetNeibFace(3,0,0,PIC::Mesh::mesh);
            }
            else if (MatrixRowNonZeroElementTable[ii].j<0) {
              MatrixRowNonZeroElementTable[ii].j+=_BLOCK_CELLS_Y_;
              MatrixRowNonZeroElementTable[ii].Node=MatrixRowNonZeroElementTable[ii].Node->GetNeibFace(2,0,0,PIC::Mesh::mesh);
            }

            if ((MatrixRowNonZeroElementTable[ii].k>=kMax) && (MatrixRowNonZeroElementTable[ii].Node->GetNeibFace(5,0,0,PIC::Mesh::mesh)!=NULL)) {
              MatrixRowNonZeroElementTable[ii].k-=_BLOCK_CELLS_Z_;
              MatrixRowNonZeroElementTable[ii].Node=MatrixRowNonZeroElementTable[ii].Node->GetNeibFace(5,0,0,PIC::Mesh::mesh);
            }
            else if (MatrixRowNonZeroElementTable[ii].k<0) {
              MatrixRowNonZeroElementTable[ii].k+=_BLOCK_CELLS_Z_;
              MatrixRowNonZeroElementTable[ii].Node=MatrixRowNonZeroElementTable[ii].Node->GetNeibFace(4,0,0,PIC::Mesh::mesh);
            }

            if (MatrixRowNonZeroElementTable[ii].Node==NULL) {
              exit(__LINE__,__FILE__,"Error: the block is not found");
            }

            cCornerNode* CornerNode=MatrixRowNonZeroElementTable[ii].Node->block->GetCornerNode(PIC::Mesh::mesh->getCornerNodeLocalNumber(MatrixRowNonZeroElementTable[ii].i,MatrixRowNonZeroElementTable[ii].j,MatrixRowNonZeroElementTable[ii].k));

            if (CornerNode==NULL) exit(__LINE__,__FILE__,"Error: something is wrong");

            el=NewRow->Elements+ii;
            el_data=NewRow->ElementDataTable+ii;

            el->CornerNode=CornerNode;
            el->CornerNodeID=PIC::Mesh::mesh->getCornerNodeLocalNumber(MatrixRowNonZeroElementTable[ii].i,MatrixRowNonZeroElementTable[ii].j,MatrixRowNonZeroElementTable[ii].k);
            el_data->UnknownVectorIndex=-1;
            el_data->MatrixElementValue=MatrixRowNonZeroElementTable[ii].MatrixElementValue;
            el->node=MatrixRowNonZeroElementTable[ii].Node;
            el_data->iVar=MatrixRowNonZeroElementTable[ii].iVar;

            el->MatrixElementParameterTableLength=MatrixRowNonZeroElementTable[ii].MatrixElementParameterTableLength;
            memcpy(el->MatrixElementParameterTable,MatrixRowNonZeroElementTable[ii].MatrixElementParameterTable,el->MatrixElementParameterTableLength*sizeof(double));

            el->MatrixElementSupportTableLength=MatrixRowNonZeroElementTable[ii].MatrixElementSupportTableLength;
            memcpy(el->MatrixElementSupportTable,MatrixRowNonZeroElementTable[ii].MatrixElementSupportTable,el->MatrixElementSupportTableLength*sizeof(void*));

            //count the number of the element that are needed
            if (CornerNode->LinearSolverUnknownVectorIndex==-1) {
              RecvDataPointCounter[MatrixRowNonZeroElementTable[ii].Node->Thread]++;
              CornerNode->LinearSolverUnknownVectorIndex=-2;
            }
          }
        }
        else {
          //the point is outside of the domain
          NewRow=MatrixRowStack.newElement();

          NewRow->i=-1000000,NewRow->j=-1000000,NewRow->k=-1000000;
          NewRow->iVar=iVar;
          NewRow->node=node;
          NewRow->Rhs=1.0;
          NewRow->CornerNode=NULL;
          NewRow->nNonZeroElements=1;
          NewRow->RhsSupportLength_CornerNodes=0;
          NewRow->RhsSupportLength_CenterNodes=0;

          NewRow->Elements[0].CornerNode=NULL;
          NewRow->Elements[0].CornerNodeID=-1;
          NewRow->ElementDataTable[0].UnknownVectorIndex=-1; //Index;
          NewRow->ElementDataTable[0].MatrixElementValue=1.0;
          NewRow->Elements[0].node=node;
          NewRow->ElementDataTable[0].iVar=iVar;

          //add the new row to the matrix
          if (MatrixRowListFirst==NULL) {
            MatrixRowListLast=NewRow;
            MatrixRowListFirst=NewRow;
            NewRow->next=NULL;
          }
          else {
            MatrixRowListLast->next=NewRow;
            MatrixRowListLast=NewRow;
            NewRow->next=NULL;
          }

        }
      }
    }
  }
}

template <class cCornerNode, int NodeUnknownVariableVectorLength,int MaxStencilLength,
int MaxRhsSupportLength_CornerNodes,int MaxRhsSupportLength_CenterNodes,
int MaxMatrixElementParameterTableLength,int MaxMatrixElementSupportTableLength>
_TARGET_HOST_
void cLinearSystemCornerNode<cCornerNode, NodeUnknownVariableVectorLength,MaxStencilLength,MaxRhsSupportLength_CornerNodes,MaxRhsSupportLength_CenterNodes,MaxMatrixElementParameterTableLength,MaxMatrixElementSupportTableLength>::BuildMatrix(void(*f)(int i,int j,int k,int iVar,cMatrixRowNonZeroElementTable* Set,int& NonZeroElementsFound,double& Rhs,cRhsSupportTable* RhsSupportTable_CornerNodes,int &RhsSupportLength_CornerNodes,cRhsSupportTable* RhsSupportTable_CenterNodes,int &RhsSupportLength_CenterNodes,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node)) {
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;
  int i,j,k,thread,nLocalNode;
  int iRow=0;
  cMatrixRow* Row;
  cStencilElement* el;
  cStencilElementData* el_data;
  cCornerNode* CornerNode;

//  int Debug=0;

  nRealBlocks=0;

  //increment the matrix modification counter 
  nMartixModifications++;

  //reset indexing of the nodes
  Reset();

  //save the user-defined build stencil function
  GetStencil=f;

  //allocate the counter of the data points to be recieve
  int *RecvDataPointCounter=new int [PIC::nTotalThreads];
  for (thread=0;thread<PIC::nTotalThreads;thread++) RecvDataPointCounter[thread]=0;

  list<cLinearSystemCornerNodeDataRequestListElement> *DataRequestList=new list<cLinearSystemCornerNodeDataRequestListElement> [PIC::nTotalThreads];

  //build the matrix
  for (nLocalNode=0;nLocalNode<PIC::DomainBlockDecomposition::nLocalBlocks;nLocalNode++) {
    node=PIC::DomainBlockDecomposition::BlockTable[nLocalNode];
    if (node->block==NULL) continue;
    //in case of the periodic boundary condition it is only the points that are inside the "real" computational domain that are considered
    if (_PIC_BC__PERIODIC_MODE_==_PIC_BC__PERIODIC_MODE_ON_) {
      bool BoundaryBlock=false;

      for (int iface=0;iface<6;iface++) if (node->GetNeibFace(iface,0,0,PIC::Mesh::mesh)==NULL) {
        //the block is at the domain boundary, and thresefor it is a 'ghost' block that is used to impose the periodic boundary conditions
        BoundaryBlock=true;
        break;
      }

      if (BoundaryBlock==true) continue;
    }

    nRealBlocks++;

    //the limits are correct: the point i==_BLOCK_CELLS_X_ belongs to the onother block
    int kMax=_BLOCK_CELLS_Z_,jMax=_BLOCK_CELLS_Y_,iMax=_BLOCK_CELLS_X_;
    int iStop,jStop,kStop;

    if ((_PIC_COUPLER_MODE_ != _PIC_COUPLER_MODE__SWMF_) && (_PIC_BC__PERIODIC_MODE_ != _PIC_BC__PERIODIC_MODE_ON_)) { 
      // Touches global domain boundary on MAX faces?
      const bool touchXmax = (node->GetNeibFace(0*2 + 1, 0, 0, PIC::Mesh::mesh) == nullptr);
      const bool touchYmax = (node->GetNeibFace(1*2 + 1, 0, 0, PIC::Mesh::mesh) == nullptr);
      const bool touchZmax = (node->GetNeibFace(2*2 + 1, 0, 0, PIC::Mesh::mesh) == nullptr);

      // Increment only if stop==max (prevents >+1 when corners belong to multiple edges).
      if (touchXmax && iStop == iMax) ++iStop;
      if (touchYmax && jStop == jMax) ++jStop;
      if (touchZmax && kStop == kMax) ++kStop;
    }
    else {
      iStop=iMax,jStop=jMax,kStop=kMax;
    }

    // Build over the (possibly-extended) domain
    for (k = 0; k < kStop; ++k)
    for (j = 0; j < jStop; ++j)
    for (i = 0; i < iStop; ++i) {
      for (int iVar=0;iVar<NodeUnknownVariableVectorLength;iVar++) {
        int NonZeroElementsFound;
        double rhs;

        //create the new Row entry
        cMatrixRow* NewRow=MatrixRowStack.newElement();

        //call the user defined function to determine the non-zero elements of the matrix
        NewRow->RhsSupportLength_CornerNodes=0;
        NewRow->RhsSupportLength_CenterNodes=0;
        rhs=0.0;

        for (int ii=0;ii<MaxStencilLength;ii++) MatrixRowNonZeroElementTable[ii].MatrixElementParameterTableLength=0,MatrixRowNonZeroElementTable[ii].MatrixElementSupportTableLength=0;

        f(i,j,k,iVar,MatrixRowNonZeroElementTable,NonZeroElementsFound,rhs,NewRow->RhsSupportTable_CornerNodes,NewRow->RhsSupportLength_CornerNodes,NewRow->RhsSupportTable_CenterNodes,NewRow->RhsSupportLength_CenterNodes,node);

        if (NonZeroElementsFound==0) {
          //the point is not included into the matrix
          MatrixRowStack.deleteElement(NewRow);
          continue;
        }

        if (NonZeroElementsFound>MaxStencilLength) {
          exit(__LINE__,__FILE__,"Error: NonZeroElementsFound>=nMaxMatrixNonzeroElement; Need to increase the value of nMaxMatrixNonzeroElement");
        }

        if (NewRow->RhsSupportLength_CenterNodes>MaxRhsSupportLength_CenterNodes) {
          exit(__LINE__,__FILE__,"Error: NewRow->RhsSupportLength_CenterNodes>MaxRhsSupportLength_CenterNodes; Need to increase the value of MaxRhsSupportLength_CenterNodes");
        }

        if (NewRow->RhsSupportLength_CornerNodes>MaxRhsSupportLength_CornerNodes) {
          exit(__LINE__,__FILE__,"Error: NewRow->RhsSupportLength_CornerNodes>MaxRhsSupportLength_CornerNodes; Need to increase the value of MaxRhsSupportLength_CornerNodes");
        }

        //scan through the found stencil and correct blocks and indexing is needed
        for (int ii=0;ii<NonZeroElementsFound;ii++) {
          MatrixRowNonZeroElementTable[ii].Node=node;

          if ((MatrixRowNonZeroElementTable[ii].i>=iMax) && (MatrixRowNonZeroElementTable[ii].Node->GetNeibFace(1,0,0,PIC::Mesh::mesh)!=NULL)) {
            MatrixRowNonZeroElementTable[ii].i-=_BLOCK_CELLS_X_;
            MatrixRowNonZeroElementTable[ii].Node=MatrixRowNonZeroElementTable[ii].Node->GetNeibFace(1,0,0,PIC::Mesh::mesh);
          }
          else if (MatrixRowNonZeroElementTable[ii].i<0) {
            MatrixRowNonZeroElementTable[ii].i+=_BLOCK_CELLS_X_;
            MatrixRowNonZeroElementTable[ii].Node=MatrixRowNonZeroElementTable[ii].Node->GetNeibFace(0,0,0,PIC::Mesh::mesh);
          }

          if ((MatrixRowNonZeroElementTable[ii].j>=jMax) && (MatrixRowNonZeroElementTable[ii].Node->GetNeibFace(3,0,0,PIC::Mesh::mesh)!=NULL)) {
            MatrixRowNonZeroElementTable[ii].j-=_BLOCK_CELLS_Y_;
            MatrixRowNonZeroElementTable[ii].Node=MatrixRowNonZeroElementTable[ii].Node->GetNeibFace(3,0,0,PIC::Mesh::mesh);
          }
          else if (MatrixRowNonZeroElementTable[ii].j<0) {
            MatrixRowNonZeroElementTable[ii].j+=_BLOCK_CELLS_Y_;
            MatrixRowNonZeroElementTable[ii].Node=MatrixRowNonZeroElementTable[ii].Node->GetNeibFace(2,0,0,PIC::Mesh::mesh);
          }

          if ((MatrixRowNonZeroElementTable[ii].k>=kMax) && (MatrixRowNonZeroElementTable[ii].Node->GetNeibFace(5,0,0,PIC::Mesh::mesh)!=NULL)) {
            MatrixRowNonZeroElementTable[ii].k-=_BLOCK_CELLS_Z_;
            MatrixRowNonZeroElementTable[ii].Node=MatrixRowNonZeroElementTable[ii].Node->GetNeibFace(5,0,0,PIC::Mesh::mesh);
          }
          else if (MatrixRowNonZeroElementTable[ii].k<0) {
            MatrixRowNonZeroElementTable[ii].k+=_BLOCK_CELLS_Z_;
            MatrixRowNonZeroElementTable[ii].Node=MatrixRowNonZeroElementTable[ii].Node->GetNeibFace(4,0,0,PIC::Mesh::mesh);
          }

          if (MatrixRowNonZeroElementTable[ii].Node==NULL) {
            exit(__LINE__,__FILE__,"Error: the block is not found");
          }

          //check whether the periodic boundary conditions are in use, and the new block is on the boundary of the domain
          if (_PIC_BC__PERIODIC_MODE_==_PIC_BC__PERIODIC_MODE_ON_) {
            bool BoundaryBlock=false;

            for (int iface=0;iface<6;iface++) if (MatrixRowNonZeroElementTable[ii].Node->GetNeibFace(iface,0,0,PIC::Mesh::mesh)==NULL) {
              //the block is at the domain boundary, and thresefor it is a 'ghost' block that is used to impose the periodic boundary conditions
              BoundaryBlock=true;
              break;
            }

            MatrixRowNonZeroElementTable[ii].BoundaryNodeFlag=BoundaryBlock;
            MatrixRowNonZeroElementTable[ii].OriginalNode=MatrixRowNonZeroElementTable[ii].Node;

            if (BoundaryBlock==true) {
              //the block is at the domain boundary -> find the point located in the "real" part of the computational domain
              MatrixRowNonZeroElementTable[ii].Node=PIC::BC::ExternalBoundary::Periodic::findCorrespondingRealBlock(MatrixRowNonZeroElementTable[ii].Node);
            }
          }

        }

        //populate the new row
        NewRow->i=i,NewRow->j=j,NewRow->k=k;
        NewRow->iVar=iVar;
        NewRow->node=node;
        NewRow->Rhs=rhs;
        NewRow->CornerNode=node->block->GetCornerNode(PIC::Mesh::mesh->getCornerNodeLocalNumber(i,j,k));
        NewRow->nNonZeroElements=NonZeroElementsFound;


        if (MatrixRowListFirst==NULL) {
          MatrixRowListLast=NewRow;
          MatrixRowListFirst=NewRow;
          NewRow->next=NULL;
        }
        else {
          MatrixRowListLast->next=NewRow;
          MatrixRowListLast=NewRow;
          NewRow->next=NULL;
        }

        //add to the row non-zero elements
        for (int iElement=0;iElement<NonZeroElementsFound;iElement++) {
          el=NewRow->Elements+iElement;
          el_data=NewRow->ElementDataTable+iElement;

          if (_PIC_BC__PERIODIC_MODE_==_PIC_BC__PERIODIC_MODE_OFF_) {
            CornerNode=MatrixRowNonZeroElementTable[iElement].Node->block->GetCornerNode(PIC::Mesh::mesh->getCornerNodeLocalNumber(MatrixRowNonZeroElementTable[iElement].i,MatrixRowNonZeroElementTable[iElement].j,MatrixRowNonZeroElementTable[iElement].k));
          }
          else  {
            if (MatrixRowNonZeroElementTable[iElement].BoundaryNodeFlag==false) {
              CornerNode=MatrixRowNonZeroElementTable[iElement].Node->block->GetCornerNode(PIC::Mesh::mesh->getCornerNodeLocalNumber(MatrixRowNonZeroElementTable[iElement].i,MatrixRowNonZeroElementTable[iElement].j,MatrixRowNonZeroElementTable[iElement].k));
            }
            else {
              if (MatrixRowNonZeroElementTable[iElement].Node->Thread!=PIC::ThisThread) {
                CornerNode=MatrixRowNonZeroElementTable[iElement].OriginalNode->block->GetCornerNode(PIC::Mesh::mesh->getCornerNodeLocalNumber(MatrixRowNonZeroElementTable[iElement].i,MatrixRowNonZeroElementTable[iElement].j,MatrixRowNonZeroElementTable[iElement].k));
              }
              else {
                CornerNode=MatrixRowNonZeroElementTable[iElement].Node->block->GetCornerNode(PIC::Mesh::mesh->getCornerNodeLocalNumber(MatrixRowNonZeroElementTable[iElement].i,MatrixRowNonZeroElementTable[iElement].j,MatrixRowNonZeroElementTable[iElement].k));
              }
            }
          }

          el->CornerNode=CornerNode;
          el->CornerNodeID=PIC::Mesh::mesh->getCornerNodeLocalNumber(MatrixRowNonZeroElementTable[iElement].i,MatrixRowNonZeroElementTable[iElement].j,MatrixRowNonZeroElementTable[iElement].k);
          el_data->UnknownVectorIndex=-1;
          el_data->MatrixElementValue=MatrixRowNonZeroElementTable[iElement].MatrixElementValue;
          el->node=MatrixRowNonZeroElementTable[iElement].Node;
          el_data->iVar=MatrixRowNonZeroElementTable[iElement].iVar;

          el->MatrixElementParameterTableLength=MatrixRowNonZeroElementTable[iElement].MatrixElementParameterTableLength;
          memcpy(el->MatrixElementParameterTable,MatrixRowNonZeroElementTable[iElement].MatrixElementParameterTable,el->MatrixElementParameterTableLength*sizeof(double));

          el->MatrixElementSupportTableLength=MatrixRowNonZeroElementTable[iElement].MatrixElementSupportTableLength;
          memcpy(el->MatrixElementSupportTable,MatrixRowNonZeroElementTable[iElement].MatrixElementSupportTable,el->MatrixElementSupportTableLength*sizeof(void*));

          //count the number of the element that are needed
          if (MatrixRowNonZeroElementTable[iElement].Node->Thread==PIC::ThisThread) {
            //the point is located at the current MPI process

            if (CornerNode->LinearSolverUnknownVectorIndex==-1) {
              RecvDataPointCounter[PIC::ThisThread]++;
              CornerNode->LinearSolverUnknownVectorIndex=-2;
            }
          }
          else {
            //the data point is on another MPI process

            if (CornerNode->LinearSolverUnknownVectorIndex==-1) {
              RecvDataPointCounter[MatrixRowNonZeroElementTable[iElement].Node->Thread]++;
              CornerNode->LinearSolverUnknownVectorIndex=-2;
            }
          }
        }


      }
    }
  }


  //add 'virtual' blocks at the right boundary of the domain
  if (_PIC_BC__PERIODIC_MODE_==_PIC_BC__PERIODIC_MODE_OFF_) {
    //ProcessRightDomainBoundary(RecvDataPointCounter,f);
  }


  //allocate the data buffers for the partial vectors and link the pointers points that are inside the subdomian
  //allocate the exchange buffer for the data the needs to be recieved from other MPI processess
  int iElement;

  int *DataExchangeTableCounter=new int [PIC::nTotalThreads];
  for (thread=0;thread<PIC::nTotalThreads;thread++) DataExchangeTableCounter[thread]=0;

  for (iRow=0,Row=MatrixRowListFirst;Row!=NULL;iRow++,Row=Row->next) {

    for (iElement=0;iElement<Row->nNonZeroElements;iElement++) {
      el=Row->Elements+iElement;
      el_data=Row->ElementDataTable+iElement;

      if (el->CornerNode!=NULL) {
        if (el->CornerNode->LinearSolverUnknownVectorIndex==-2) {
          //the node is still not linked to the 'partial data vector'
          el->CornerNode->LinearSolverUnknownVectorIndex=DataExchangeTableCounter[el->node->Thread]++;

          if (el->node->Thread!=PIC::ThisThread) {
            //add the point to the data request list
            cLinearSystemCornerNodeDataRequestListElement DataRequestListElement;
            cAMRnodeID NodeID;

            PIC::Mesh::mesh->GetAMRnodeID(NodeID,el->node);
            DataRequestListElement.CornerNodeID=el->CornerNodeID;
            DataRequestListElement.NodeID=NodeID;

            DataRequestList[el->node->Thread].push_back(DataRequestListElement);
          }
        }

        el_data->UnknownVectorIndex=el->CornerNode->LinearSolverUnknownVectorIndex;
        el_data->Thread=el->node->Thread;
nMartixModifications++;
      }
      else {
        //the point is within the 'virtual' block
        el_data->UnknownVectorIndex=DataExchangeTableCounter[el->node->Thread]++;
        el_data->Thread=el->node->Thread;
nMartixModifications++;
      }
    }
  }

  //re-count the 'internal' corner nodes so they follow the i,j,k order as well as rows
  for (Row=MatrixRowListFirst;Row!=NULL;Row=Row->next) if (Row->CornerNode!=NULL) Row->CornerNode->LinearSolverUnknownVectorIndex=-1;

  for (iRow=0,Row=MatrixRowListFirst;Row!=NULL;Row=Row->next) {
    if (Row->CornerNode!=NULL) {
      if (Row->iVar==0) {
        if (Row->CornerNode->LinearSolverUnknownVectorIndex>=0) exit(__LINE__,__FILE__,"Error: a courner node is double counted");
        Row->CornerNode->LinearSolverUnknownVectorIndex=iRow;
      }
    }
    else Row->ElementDataTable[0].UnknownVectorIndex=iRow; //rows corresponding to the virtual blocks has only one non-zero element

    if (Row->iVar==NodeUnknownVariableVectorLength-1) iRow++; //all rows pointing to the same 'CornerNode' have the same 'iRow'
  }

  //verify that all all corner nodes have been counted
  for (Row=MatrixRowListFirst;Row!=NULL;iRow++,Row=Row->next) for (iElement=0;iElement<Row->nNonZeroElements;iElement++) {
    int n;

    if (Row->Elements[iElement].CornerNode!=NULL) {
      if ((n=Row->Elements[iElement].CornerNode->LinearSolverUnknownVectorIndex)<0) {
        exit(__LINE__,__FILE__,"Error: an uncounted corner node has been found");
      }

      Row->ElementDataTable[iElement].UnknownVectorIndex=n;
    }
  }


  //exchange the request lists
  int From,To;
  int *nGlobalDataPointTable=new int [PIC::nTotalThreads*PIC::nTotalThreads];
  cLinearSystemCornerNodeDataRequestListElement* ExchangeList;

  MPI_Allgather(DataExchangeTableCounter,PIC::nTotalThreads,MPI_INT,nGlobalDataPointTable,PIC::nTotalThreads,MPI_INT,MPI_GLOBAL_COMMUNICATOR);

  SendExchangeBufferLength=new int [PIC::nTotalThreads];
  SendExchangeBufferElementIndex=new int* [PIC::nTotalThreads];

  RecvExchangeBufferLength=new int [PIC::nTotalThreads];

  for (thread=0;thread<PIC::nTotalThreads;thread++) {
    SendExchangeBufferLength[thread]=0,SendExchangeBufferElementIndex[thread]=0;
    RecvExchangeBufferLength[thread]=0;
  }


  //create the lists
  for (To=0;To<PIC::nTotalThreads;To++) for (From=0;From<PIC::nTotalThreads;From++) if ( (To!=From) && (nGlobalDataPointTable[From+To*PIC::nTotalThreads]!=0) && ((To==PIC::ThisThread)||(From==PIC::ThisThread)) ) {
    ExchangeList=new cLinearSystemCornerNodeDataRequestListElement [nGlobalDataPointTable[From+To*PIC::nTotalThreads]];

    if (PIC::ThisThread==To) {
      list<cLinearSystemCornerNodeDataRequestListElement>::iterator ptr;

      for (i=0,ptr=DataRequestList[From].begin();ptr!=DataRequestList[From].end();ptr++,i++) ExchangeList[i]=*ptr;
      MPI_Send(ExchangeList,nGlobalDataPointTable[From+To*PIC::nTotalThreads]*sizeof(cLinearSystemCornerNodeDataRequestListElement),MPI_CHAR,From,0,MPI_GLOBAL_COMMUNICATOR);

      //create the recv list
      RecvExchangeBufferLength[From]=nGlobalDataPointTable[From+To*PIC::nTotalThreads];
    }
    else {
      MPI_Status status;

      MPI_Recv(ExchangeList,nGlobalDataPointTable[From+To*PIC::nTotalThreads]*sizeof(cLinearSystemCornerNodeDataRequestListElement),MPI_CHAR,To,0,MPI_GLOBAL_COMMUNICATOR,&status);

      //unpack the SendDataList
      SendExchangeBufferLength[To]=nGlobalDataPointTable[From+To*PIC::nTotalThreads];
      SendExchangeBufferElementIndex[To]=new int[SendExchangeBufferLength[To]];

      for (int ii=0;ii<SendExchangeBufferLength[To];ii++) {
        node=PIC::Mesh::mesh->findAMRnodeWithID(ExchangeList[ii].NodeID);
        CornerNode=node->block->GetCornerNode(ExchangeList[ii].CornerNodeID);

        SendExchangeBufferElementIndex[To][ii]=CornerNode->LinearSolverUnknownVectorIndex;
      }
    }

    delete [] ExchangeList;
  }

  //deallocate 'nGlobalDataPointTable'
  delete [] nGlobalDataPointTable;
  delete [] DataExchangeTableCounter;

  delete [] DataRequestList;
  delete [] RecvDataPointCounter;

  //cleate a RowTable;
  for (MatrixRowTableLength=0,Row=MatrixRowListFirst;Row!=NULL;Row=Row->next) MatrixRowTableLength++;

  if (MatrixRowTable!=NULL) {
    //delete [] MatrixRowTable;
    amps_free_managed(MatrixRowTable);
    MatrixRowTable=NULL;
  }

  //MatrixRowTable=new cMatrixRow *[MatrixRowTableLength];
  amps_new_managed(MatrixRowTable,MatrixRowTableLength);

  for (MatrixRowTableLength=0,Row=MatrixRowListFirst;Row!=NULL;Row=Row->next) MatrixRowTable[MatrixRowTableLength++]=Row;
}

template <class cCornerNode, int NodeUnknownVariableVectorLength,int MaxStencilLength,
int MaxRhsSupportLength_CornerNodes,int MaxRhsSupportLength_CenterNodes,
int MaxMatrixElementParameterTableLength,int MaxMatrixElementSupportTableLength>
void cLinearSystemCornerNode<cCornerNode, NodeUnknownVariableVectorLength,MaxStencilLength,MaxRhsSupportLength_CornerNodes,MaxRhsSupportLength_CenterNodes,MaxMatrixElementParameterTableLength,MaxMatrixElementSupportTableLength>::ExchangeIntermediateUnknownsData(double *x,double** &LocalRecvExchangeBufferTable) {
  int To,From;

//  LocalRecvExchangeBufferTable=new double *[PIC::nTotalThreads];
  
  MPI_Datatype *unknown_vector_type_table=new MPI_Datatype[PIC::nTotalThreads];
  
  if (LocalRecvExchangeBufferTable==NULL) {
    //LocalRecvExchangeBufferTable=new double *[PIC::nTotalThreads]; 
    amps_new_managed(LocalRecvExchangeBufferTable,PIC::nTotalThreads);

  for (int thread=0;thread<PIC::nTotalThreads;thread++) {
   LocalRecvExchangeBufferTable[thread]=NULL;

    if (NodeUnknownVariableVectorLength*RecvExchangeBufferLength[thread]>0) {
      LocalRecvExchangeBufferTable[thread]=new double [NodeUnknownVariableVectorLength*RecvExchangeBufferLength[thread]];
    }
  } 
}


  MPI_Request RecvRequestTable[PIC::nTotalThreads],SendRequestTable[PIC::nTotalThreads];
  int RecvRequestTableLength=0,SendRequestTableLength=0;

  //Initiate Recv operations
  for (From=0;From<PIC::nTotalThreads;From++) if ((From!=PIC::ThisThread)&&(RecvExchangeBufferLength[From]!=0)) {
    MPI_Irecv(LocalRecvExchangeBufferTable[From],NodeUnknownVariableVectorLength*RecvExchangeBufferLength[From],MPI_DOUBLE,From,0,MPI_GLOBAL_COMMUNICATOR,RecvRequestTable+RecvRequestTableLength);
    RecvRequestTableLength++;
  }

  //Init Send operations
  for (To=0;To<PIC::nTotalThreads;To++) if ((To!=PIC::ThisThread)&&(SendExchangeBufferLength[To]!=0)) {
    //this if the sending thread
    int i,offset=0;
    int Length=SendExchangeBufferLength[To];
    int *IndexTable=SendExchangeBufferElementIndex[To];

    int length_table[Length];
    int offset_table[Length];

    for (i=0;i<Length;i++) {
      length_table[i]=NodeUnknownVariableVectorLength;
      offset_table[i]=NodeUnknownVariableVectorLength*IndexTable[i]; 
    }

    MPI_Type_indexed(Length, length_table, offset_table, MPI_DOUBLE, unknown_vector_type_table+To);
    MPI_Type_commit(unknown_vector_type_table+To);

    for (i=0;i<Length;i++) {
      offset+=NodeUnknownVariableVectorLength;
    }

    if (offset!=NodeUnknownVariableVectorLength*SendExchangeBufferLength[To]) exit(__LINE__,__FILE__,"Error: out of limit");

    MPI_Isend(x,1,unknown_vector_type_table[To],To,0,MPI_GLOBAL_COMMUNICATOR,SendRequestTable+SendRequestTableLength);

    MPI_Type_free(unknown_vector_type_table+To); 
    SendRequestTableLength++;
  }

  //wait for completing the send/rect operations
  MPI_Waitall(RecvRequestTableLength,RecvRequestTable,MPI_STATUSES_IGNORE);
  MPI_Waitall(SendRequestTableLength,SendRequestTable,MPI_STATUSES_IGNORE);

  //delete temporary buffers
  delete [] unknown_vector_type_table;
}

template <class cCornerNode, int NodeUnknownVariableVectorLength,int MaxStencilLength,
int MaxRhsSupportLength_CornerNodes,int MaxRhsSupportLength_CenterNodes,
int MaxMatrixElementParameterTableLength,int MaxMatrixElementSupportTableLength>
_TARGET_HOST_ _TARGET_DEVICE_
void cLinearSystemCornerNode<cCornerNode, NodeUnknownVariableVectorLength,MaxStencilLength,MaxRhsSupportLength_CornerNodes,MaxRhsSupportLength_CenterNodes,MaxMatrixElementParameterTableLength,MaxMatrixElementSupportTableLength>::DeleteDataBuffers() {
  #ifdef __CUDA_ARCH__ 
  using namespace PIC::GPU; 
  #else 
  using namespace PIC::CPU; 
  #endif

  if (RecvExchangeBufferLength!=NULL) {
    //clear the row stack
    MatrixRowStack.resetStack();

    //deallocate the allocated data buffers
    for (int thread=0;thread<nTotalThreads;thread++) {

      if (SendExchangeBufferElementIndex!=NULL) if (SendExchangeBufferElementIndex[thread]!=NULL) {
        delete [] SendExchangeBufferElementIndex[thread];
      }
    }

    if (RecvExchangeBufferLength!=NULL) {
      delete [] RecvExchangeBufferLength;
    }

    RecvExchangeBufferLength=NULL;

    if (SendExchangeBufferElementIndex!=NULL) {
      delete [] SendExchangeBufferLength;
      delete [] SendExchangeBufferElementIndex;
    }

    SendExchangeBufferLength=NULL,SendExchangeBufferElementIndex=NULL;

    if (SubdomainPartialRHS!=NULL) {
      delete [] SubdomainPartialRHS;
      delete [] SubdomainPartialUnknownsVector;
    }

    MatrixRowListFirst=NULL,MatrixRowListLast=NULL;
    SubdomainPartialRHS=NULL,SubdomainPartialUnknownsVector=NULL;
  }
}

template <class cCornerNode, int NodeUnknownVariableVectorLength,int MaxStencilLength,
int MaxRhsSupportLength_CornerNodes,int MaxRhsSupportLength_CenterNodes,
int MaxMatrixElementParameterTableLength,int MaxMatrixElementSupportTableLength>
_TARGET_HOST_ _TARGET_DEVICE_
void cLinearSystemCornerNode<cCornerNode, NodeUnknownVariableVectorLength,MaxStencilLength,MaxRhsSupportLength_CornerNodes,MaxRhsSupportLength_CenterNodes,MaxMatrixElementParameterTableLength,MaxMatrixElementSupportTableLength>::Reset() {

  cAmpsMesh<PIC::Mesh::cDataCornerNode,PIC::Mesh::cDataCenterNode,PIC::Mesh::cDataBlockAMR>  *mesh=PIC::Mesh::mesh;


  Reset(mesh->rootTree);
}

template <class cCornerNode, int NodeUnknownVariableVectorLength,int MaxStencilLength,
int MaxRhsSupportLength_CornerNodes,int MaxRhsSupportLength_CenterNodes,
int MaxMatrixElementParameterTableLength,int MaxMatrixElementSupportTableLength>
_TARGET_HOST_ _TARGET_DEVICE_
void cLinearSystemCornerNode<cCornerNode, NodeUnknownVariableVectorLength,MaxStencilLength,MaxRhsSupportLength_CornerNodes,MaxRhsSupportLength_CenterNodes,MaxMatrixElementParameterTableLength,MaxMatrixElementSupportTableLength>::Reset(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {


  cAmpsMesh<PIC::Mesh::cDataCornerNode,PIC::Mesh::cDataCenterNode,PIC::Mesh::cDataBlockAMR>  *mesh=PIC::Mesh::mesh;


  if ((startNode==mesh->rootTree)&&(RecvExchangeBufferLength!=NULL)) DeleteDataBuffers();

  if (startNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
    PIC::Mesh::cDataBlockAMR *block=NULL;
    PIC::Mesh::cDataCornerNode *CornerNode=NULL;

    if ((block=startNode->block)!=NULL) {
      for (int i=0;i<_BLOCK_CELLS_X_;i++) for (int j=0;j<_BLOCK_CELLS_Y_;j++) for (int k=0;k<_BLOCK_CELLS_Z_;k++) {
        if ((CornerNode=block->GetCornerNode(mesh->getCornerNodeLocalNumber(i,j,k)))!=NULL) {
          CornerNode->LinearSolverUnknownVectorIndex=-1;
        }
      }
    }

  }
  else for (int i=0;i<(1<<DIM);i++) if (startNode->downNode[i]!=NULL) Reset(startNode->downNode[i]);
}


//update non-zero coefficients of the matrix
template <class cCornerNode, int NodeUnknownVariableVectorLength,int MaxStencilLength,
int MaxRhsSupportLength_CornerNodes,int MaxRhsSupportLength_CenterNodes,
int MaxMatrixElementParameterTableLength,int MaxMatrixElementSupportTableLength>
void cLinearSystemCornerNode<cCornerNode, NodeUnknownVariableVectorLength,MaxStencilLength,MaxRhsSupportLength_CornerNodes,MaxRhsSupportLength_CenterNodes,MaxMatrixElementParameterTableLength,MaxMatrixElementSupportTableLength>::UpdateMatrixNonZeroCoefficients(void (*UpdateMatrixRow)(cMatrixRow*)) {

  #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
  #pragma omp parallel for schedule(dynamic,1)
  #endif
  for (int irow=0;irow<MatrixRowTableLength;irow++) {
    UpdateMatrixRow(MatrixRowTable[irow]);
  }
}


template <class cCornerNode, int NodeUnknownVariableVectorLength,int MaxStencilLength,
int MaxRhsSupportLength_CornerNodes,int MaxRhsSupportLength_CenterNodes,
int MaxMatrixElementParameterTableLength,int MaxMatrixElementSupportTableLength>
void cLinearSystemCornerNode<cCornerNode, NodeUnknownVariableVectorLength,MaxStencilLength,MaxRhsSupportLength_CornerNodes,MaxRhsSupportLength_CenterNodes,MaxMatrixElementParameterTableLength,MaxMatrixElementSupportTableLength>::UpdateRhs(double (*fSetRhs)(int,cRhsSupportTable*,int,cRhsSupportTable*,int)) {

  #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
  #pragma omp parallel for schedule(dynamic,1) default(none) shared(fSetRhs)
  #endif
  for (int irow=0;irow<MatrixRowTableLength;irow++) {
    cMatrixRow* row=MatrixRowTable[irow];

    if ((row->RhsSupportLength_CornerNodes!=0)||(row->RhsSupportLength_CenterNodes!=0)) {
      row->Rhs=fSetRhs(row->iVar,row->RhsSupportTable_CornerNodes,row->RhsSupportLength_CornerNodes,row->RhsSupportTable_CenterNodes,row->RhsSupportLength_CenterNodes);
    }
  }
}




//CUDA version of the vector multiplication function 
#if _CUDA_MODE_ == _ON_ 
template <class cCornerNode, int NodeUnknownVariableVectorLength,int MaxStencilLength,
int MaxRhsSupportLength_CornerNodes,int MaxRhsSupportLength_CenterNodes,
int MaxMatrixElementParameterTableLength,int MaxMatrixElementSupportTableLength>
_TARGET_HOST_
void cLinearSystemCornerNode<cCornerNode, NodeUnknownVariableVectorLength,MaxStencilLength,MaxRhsSupportLength_CornerNodes,MaxRhsSupportLength_CenterNodes,MaxMatrixElementParameterTableLength,MaxMatrixElementSupportTableLength>::MultiplyVector(double *p,double *x,int length) {


  static double **RecvExchangeBufferTable=NULL;
  static int nMartixModificationsLocal=-1;  

  static double *x_local=NULL;
  static int length_local=-1;

  if (x_local==NULL) {
    //x_local=new double [length];
    amps_new_managed(x_local,length);
  }
  else if (length_local!=length) {
    //delete [] x_local;
    //x_local=new double [length];
    
    amps_free_managed(x_local);
    amps_new_managed(x_local,length); 
  } 

  memcpy(x_local,x,length*sizeof(double));
  

  if ((nMartixModificationsLocal!=nMartixModifications)&&(RecvExchangeBufferTable!=NULL)) {
    RecvExchangeBufferTable[PIC::ThisThread]=NULL;

    for (int thread=0;thread<PIC::nTotalThreads;thread++) if (RecvExchangeBufferTable[thread]!=NULL) {
     // delete [] RecvExchangeBufferTable[thread];
     amps_free_managed(RecvExchangeBufferTable[thread]);
    }

    //delete [] RecvExchangeBufferTable;
    amps_free_managed(RecvExchangeBufferTable); 

    RecvExchangeBufferTable=NULL;
  }

  ExchangeIntermediateUnknownsData(x,RecvExchangeBufferTable);
  RecvExchangeBufferTable[PIC::ThisThread]=x_local;


  auto UpdateRhsPtr = [] _TARGET_DEVICE_ (int MatrixRowTableLength, cMatrixRow** MatrixRowTable,double** RecvExchangeBufferTable) {
    int id=blockIdx.x*blockDim.x+threadIdx.x;
    int increment=gridDim.x*blockDim.x;

    for (int irow=id;irow<MatrixRowTableLength;irow+=increment) {
      int iElementMax=MatrixRowTable[irow]->nNonZeroElements;
      cStencilElementData *ElementDataTable=MatrixRowTable[irow]->ElementDataTable;

      for (int iElement=0;iElement<iElementMax;iElement++) {
        cStencilElementData *data=ElementDataTable+iElement;

        data->rhs=&RecvExchangeBufferTable[data->Thread][data->iVar+NodeUnknownVariableVectorLength*data->UnknownVectorIndex];
      }
    }
  }; 


  if ((nMartixModificationsLocal!=nMartixModifications)||(length_local!=length)) {
    //UpdateRhsPtr(MatrixRowTableLength,MatrixRowTable,RecvExchangeBufferTable);
    kernel_3<<<_CUDA_BLOCKS_,_CUDA_THREADS_>>>(UpdateRhsPtr,MatrixRowTableLength,MatrixRowTable,RecvExchangeBufferTable);
    cudaDeviceSynchronize();

    nMartixModificationsLocal=nMartixModifications;
    length_local=length; 
  }


  auto ProcessAllRows = [] _TARGET_DEVICE_ (double* p,int MatrixRowTableLength,cMatrixRow** MatrixRowTable,int length) { 
    int id=blockIdx.x*blockDim.x+threadIdx.x;
    int increment=gridDim.x*blockDim.x;

    for (int irow=id;irow<MatrixRowTableLength;irow+=increment) {
      int iElementMax=MatrixRowTable[irow]->nNonZeroElements;
      cStencilElementData *ElementDataTable=MatrixRowTable[irow]->ElementDataTable;
    
      double res=0.0;
      int iElement=0;
      cStencilElementData *data,*data_next,*data_next_next=NULL;
      double *u_vect,*u_vect_next;

      for (iElement=0;iElement<iElementMax;iElement++) {
        data=ElementDataTable+iElement;
        res+=data->MatrixElementValue*(*data->rhs);
      }

       #if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
       if (irow>=length) exit(__LINE__,__FILE__,"Error: out of bound");
       #endif

       p[irow]=res;
      __syncthreads();
    } 
  }; 

  double *p_gpu=NULL;

  cudaMalloc(&p_gpu,MatrixRowTableLength*sizeof(double)); 

  kernel_4<<<_CUDA_BLOCKS_,_CUDA_THREADS_>>>(ProcessAllRows,p_gpu,MatrixRowTableLength,MatrixRowTable,length); 
  cudaDeviceSynchronize();

  cudaMemcpy(p,p_gpu,MatrixRowTableLength*sizeof(double),cudaMemcpyDeviceToHost);
  cudaFree(p_gpu);
}

//CPU version of the vector multipllication function
#else //__CUDA_ARCH__
#if _PIC_MATMUL_MPI_MULTITHREAD_ == _PIC_MODE_ON_
template <class cCornerNode, int NodeUnknownVariableVectorLength,int MaxStencilLength,
int MaxRhsSupportLength_CornerNodes,int MaxRhsSupportLength_CenterNodes,
int MaxMatrixElementParameterTableLength,int MaxMatrixElementSupportTableLength>
_TARGET_DEVICE_
void cLinearSystemCornerNode<cCornerNode, NodeUnknownVariableVectorLength,MaxStencilLength,MaxRhsSupportLength_CornerNodes,MaxRhsSupportLength_CenterNodes,MaxMatrixElementParameterTableLength,MaxMatrixElementSupportTableLength>::MultiplyVector(double *p,double *x,int length) {

  static double **RecvExchangeBufferTable=NULL;
  static int nMartixModificationsLocal=-1;

  static double *x_local=NULL;
  static int length_local=-1;

  if (x_local==NULL) {
    x_local=new double [length];
  }
  else if (length_local!=length) {
    delete [] x_local;

    x_local=new double [length];
  }

  memcpy(x_local,x,length*sizeof(double));


  if ((nMartixModificationsLocal!=nMartixModifications)&&(RecvExchangeBufferTable!=NULL)) {
    RecvExchangeBufferTable[PIC::ThisThread]=NULL;

    for (int thread=0;thread<PIC::nTotalThreads;thread++) if (RecvExchangeBufferTable[thread]!=NULL) {
      delete [] RecvExchangeBufferTable[thread];
    }

    delete [] RecvExchangeBufferTable;
    RecvExchangeBufferTable=NULL;
  }

  ExchangeIntermediateUnknownsData(x,RecvExchangeBufferTable);
  RecvExchangeBufferTable[PIC::ThisThread]=x_local;

  if ((nMartixModificationsLocal!=nMartixModifications)||(length_local!=length)) {
    for (int irow=0;irow<MatrixRowTableLength;irow++) {
      int iElementMax=MatrixRowTable[irow]->nNonZeroElements;
      cStencilElementData *ElementDataTable=MatrixRowTable[irow]->ElementDataTable;

      for (int iElement=0;iElement<iElementMax;iElement++) {
        cStencilElementData *data=ElementDataTable+iElement;

        data->rhs=&RecvExchangeBufferTable[data->Thread][data->iVar+NodeUnknownVariableVectorLength*data->UnknownVectorIndex];
      }
    }

    nMartixModificationsLocal=nMartixModifications;
    length_local=length;
  }


  auto ProcessRow = [this,length,RecvExchangeBufferTable,p]  (int thread_id_table_size) {
    static Thread::Sync::cBarrier barrier(thread_id_table_size);

    double **LocalRecvExchangeBufferTable;

    LocalRecvExchangeBufferTable=new double*[PIC::nTotalThreads];
    for (int thread=0;thread<PIC::nTotalThreads;thread++) LocalRecvExchangeBufferTable[thread]=RecvExchangeBufferTable[thread];

    int irow,increment,irow_max_thread,iset=0;

    if (increment==0) increment=MatrixRowTableLength/(5*thread_id_table_size);
    if (increment==0) increment=MatrixRowTableLength/thread_id_table_size;
    if (increment==0) increment=1;

    atomic<int> irow_max;

    irow_max=0;

    barrier.Sync();

    do {
      irow=irow_max.fetch_add(increment);

      irow_max_thread=irow+increment;
      if (irow_max_thread>MatrixRowTableLength) irow_max_thread=MatrixRowTableLength;

      for (;irow<irow_max_thread;irow++) {

        int iElementMax=MatrixRowTable[irow]->nNonZeroElements;
        cStencilElementData *ElementDataTable=MatrixRowTable[irow]->ElementDataTable;

        double res=0.0;
        int iElement=0;
        cStencilElementData *data,*data_next,*data_next_next=NULL;
        double *u_vect,*u_vect_next;
#if _AVX_MATMUL_ == _ON_
#if _AVX_INSTRUCTIONS_USAGE_MODE_ == _AVX_INSTRUCTIONS_USAGE_MODE__256_ 
        //add most of the vector
        for (;iElement+3<iElementMax;iElement+=4) {
          alignas(32) double a[4],b[4],*r;
          __m256d av,bv,cv,rv;

          data=ElementDataTable+iElement;
          av=_mm256_set_pd((data+3)->MatrixElementValue,(data+2)->MatrixElementValue,(data+1)->MatrixElementValue,data->MatrixElementValue);

          bv=_mm256_set_pd(
              LocalRecvExchangeBufferTable[(data+3)->Thread][(data+3)->iVar+NodeUnknownVariableVectorLength*(data+3)->UnknownVectorIndex],
              LocalRecvExchangeBufferTable[(data+2)->Thread][(data+2)->iVar+NodeUnknownVariableVectorLength*(data+2)->UnknownVectorIndex],
              LocalRecvExchangeBufferTable[(data+1)->Thread][(data+1)->iVar+NodeUnknownVariableVectorLength*(data+1)->UnknownVectorIndex],
              LocalRecvExchangeBufferTable[data->Thread][data->iVar+NodeUnknownVariableVectorLength*data->UnknownVectorIndex]
          );


          cv=_mm256_mul_pd(av,bv);
          rv=_mm256_hadd_pd(cv,cv);

          r=(double*)&rv;
          res+=r[1]+r[2];
        }

#elif _AVX_INSTRUCTIONS_USAGE_MODE_ == _AVX_INSTRUCTIONS_USAGE_MODE__512_
        //add most of the vector
        for (;iElement+7<iElementMax;iElement+=8) {
          __m512d av,bv,cv;

          data=ElementDataTable+iElement;

          av=_mm512_set_pd(
              (data+7)->MatrixElementValue,(data+6)->MatrixElementValue,(data+5)->MatrixElementValue,(data+4)->MatrixElementValue,
              (data+3)->MatrixElementValue,(data+2)->MatrixElementValue,(data+1)->MatrixElementValue,data->MatrixElementValue);

          bv=_mm512_set_pd(
              LocalRecvExchangeBufferTable[(data+7)->Thread][(data+7)->iVar+NodeUnknownVariableVectorLength*(data+7)->UnknownVectorIndex],
              LocalRecvExchangeBufferTable[(data+6)->Thread][(data+6)->iVar+NodeUnknownVariableVectorLength*(data+6)->UnknownVectorIndex],
              LocalRecvExchangeBufferTable[(data+5)->Thread][(data+5)->iVar+NodeUnknownVariableVectorLength*(data+5)->UnknownVectorIndex],
              LocalRecvExchangeBufferTable[(data+4)->Thread][(data+4)->iVar+NodeUnknownVariableVectorLength*(data+4)->UnknownVectorIndex],

              LocalRecvExchangeBufferTable[(data+3)->Thread][(data+3)->iVar+NodeUnknownVariableVectorLength*(data+3)->UnknownVectorIndex],
              LocalRecvExchangeBufferTable[(data+2)->Thread][(data+2)->iVar+NodeUnknownVariableVectorLength*(data+2)->UnknownVectorIndex],
              LocalRecvExchangeBufferTable[(data+1)->Thread][(data+1)->iVar+NodeUnknownVariableVectorLength*(data+1)->UnknownVectorIndex],
              LocalRecvExchangeBufferTable[data->Thread][data->iVar+NodeUnknownVariableVectorLength*data->UnknownVectorIndex]
          );

          cv=_mm512_mul_pd(av,bv);
          res+=_mm512_reduce_add_pd(cv);
        }

        //add the rest of the vector
        for (;iElement+3<iElementMax;iElement+=4) {
          alignas(32) double a[4],b[4],*r;
          __m256d av,bv,cv,rv;

          data=ElementDataTable+iElement;

          av=_mm256_set_pd((data+3)->MatrixElementValue,(data+2)->MatrixElementValue,(data+1)->MatrixElementValue,data->MatrixElementValue);

          bv=_mm256_set_pd(
              LocalRecvExchangeBufferTable[(data+3)->Thread][(data+3)->iVar+NodeUnknownVariableVectorLength*(data+3)->UnknownVectorIndex],
              LocalRecvExchangeBufferTable[(data+2)->Thread][(data+2)->iVar+NodeUnknownVariableVectorLength*(data+2)->UnknownVectorIndex],
              LocalRecvExchangeBufferTable[(data+1)->Thread][(data+1)->iVar+NodeUnknownVariableVectorLength*(data+1)->UnknownVectorIndex],
              LocalRecvExchangeBufferTable[data->Thread][data->iVar+NodeUnknownVariableVectorLength*data->UnknownVectorIndex]
          );

          cv=_mm256_mul_pd(av,bv);
          rv=_mm256_hadd_pd(cv,cv);

          r=(double*)&rv;
          res+=r[1]+r[2];
        }
#endif
#endif
        //add the rest of the vector
        for (;iElement<iElementMax;iElement++) {
          data=ElementDataTable+iElement;
          res+=data->MatrixElementValue*LocalRecvExchangeBufferTable[data->Thread][data->iVar+NodeUnknownVariableVectorLength*data->UnknownVectorIndex];
        }

#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
        if (irow>=length) exit(__LINE__,__FILE__,"Error: out of bound");
#endif

        p[irow]=res;
      }

    }
    while  (irow_max_thread<MatrixRowTableLength);

    delete [] LocalRecvExchangeBufferTable;
  };


  //start threads
  int thread_id_table_size=4;
  std::thread tTable[thread_id_table_size];

  for (int i=1;i<thread_id_table_size;i++) {
    tTable[i]=std::thread(ProcessRow,thread_id_table_size);
  }

  ProcessRow(thread_id_table_size);

  for (int i=1;i<thread_id_table_size;i++) {
    tTable[i].join();
  }

}


#else //_PIC_MATMUL_MPI_MULTITHREAD_
template <class cCornerNode, int NodeUnknownVariableVectorLength,int MaxStencilLength,
int MaxRhsSupportLength_CornerNodes,int MaxRhsSupportLength_CenterNodes,
int MaxMatrixElementParameterTableLength,int MaxMatrixElementSupportTableLength>
_TARGET_HOST_
void cLinearSystemCornerNode<cCornerNode, NodeUnknownVariableVectorLength,MaxStencilLength,MaxRhsSupportLength_CornerNodes,MaxRhsSupportLength_CenterNodes,MaxMatrixElementParameterTableLength,MaxMatrixElementSupportTableLength>::MultiplyVector(double *p,double *x,int length) {


  static double **RecvExchangeBufferTable=NULL;
  static int nMartixModificationsLocal=-1;  

  static double *x_local=NULL;
  static int length_local=-1;

  if (x_local==NULL) {
    x_local=new double [length];
  }
  else if (length_local!=length) {
    delete [] x_local;

    x_local=new double [length];
  } 

  memcpy(x_local,x,length*sizeof(double));
  

  if ((nMartixModificationsLocal!=nMartixModifications)&&(RecvExchangeBufferTable!=NULL)) {
    RecvExchangeBufferTable[PIC::ThisThread]=NULL;

    for (int thread=0;thread<PIC::nTotalThreads;thread++) if (RecvExchangeBufferTable[thread]!=NULL) {
      delete [] RecvExchangeBufferTable[thread];
    }

    delete [] RecvExchangeBufferTable;
    RecvExchangeBufferTable=NULL;
  }

  ExchangeIntermediateUnknownsData(x,RecvExchangeBufferTable);
  RecvExchangeBufferTable[PIC::ThisThread]=x_local;

  if ((nMartixModificationsLocal!=nMartixModifications)||(length_local!=length)) {
    for (int irow=0;irow<MatrixRowTableLength;irow++) {
      int iElementMax=MatrixRowTable[irow]->nNonZeroElements;
      cStencilElementData *ElementDataTable=MatrixRowTable[irow]->ElementDataTable;

      for (int iElement=0;iElement<iElementMax;iElement++) {
        cStencilElementData *data=ElementDataTable+iElement;

        data->rhs=&RecvExchangeBufferTable[data->Thread][data->iVar+NodeUnknownVariableVectorLength*data->UnknownVectorIndex]; 
      }
    }

    nMartixModificationsLocal=nMartixModifications;
    length_local=length; 
  }


#if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
#pragma omp parallel default(none) shared (RecvExchangeBufferTable,PIC::ThisThread,PIC::nTotalThreads,p) firstprivate(length)
#endif
{


#if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
#pragma omp for schedule(dynamic,4) 
#endif
  for (int irow=0;irow<MatrixRowTableLength;irow++) {

    int iElementMax=MatrixRowTable[irow]->nNonZeroElements;
    cStencilElementData *ElementDataTable=MatrixRowTable[irow]->ElementDataTable;
    
    double res=0.0;
    int iElement=0;
    cStencilElementData *data,*data_next,*data_next_next=NULL;
    double *u_vect,*u_vect_next;

    #if _AVX_MATMUL_ == _ON_
    #if _AVX_INSTRUCTIONS_USAGE_MODE_ == _AVX_INSTRUCTIONS_USAGE_MODE__256_
    //add most of the vector
    for (;iElement+3<iElementMax;iElement+=4) {
      __m256d cv,rv;

      data=ElementDataTable+iElement;

      cv=_mm256_mul_pd(
        _mm256_set_pd((data+3)->MatrixElementValue,(data+2)->MatrixElementValue,(data+1)->MatrixElementValue,data->MatrixElementValue),         
        _mm256_set_pd(*((data+3)->rhs),*((data+2)->rhs),*((data+1)->rhs),*data->rhs)); 

      rv=_mm256_hadd_pd(cv,cv);

      res+=rv[1]+rv[2];
    }

    #elif _AVX_INSTRUCTIONS_USAGE_MODE_ == _AVX_INSTRUCTIONS_USAGE_MODE__512_
    //add most of the vector
    for (;iElement+7<iElementMax;iElement+=8) {
      __m512d av,bv,cv;

      data=ElementDataTable+iElement;

      av=_mm512_set_pd(
		      (data+7)->MatrixElementValue,(data+6)->MatrixElementValue,(data+5)->MatrixElementValue,(data+4)->MatrixElementValue,
		      (data+3)->MatrixElementValue,(data+2)->MatrixElementValue,(data+1)->MatrixElementValue,data->MatrixElementValue);  

      bv=_mm512_set_pd(
		      *((data+7)->rhs),*((data+6)->rhs),*((data+5)->rhs),*((data+4)->rhs), 
                      *((data+3)->rhs),*((data+2)->rhs),*((data+1)->rhs),*((data+0)->rhs)); 

      cv=_mm512_mul_pd(av,bv);
      res+=_mm512_reduce_add_pd(cv);
    }

    //add the rest of the vector
    for (;iElement+3<iElementMax;iElement+=4) {
      __m256d cv,rv;

      data=ElementDataTable+iElement;

      cv=_mm256_mul_pd(
        _mm256_set_pd((data+3)->MatrixElementValue,(data+2)->MatrixElementValue,(data+1)->MatrixElementValue,data->MatrixElementValue),
        _mm256_set_pd(*((data+3)->rhs),*((data+2)->rhs),*((data+1)->rhs),*data->rhs));

      rv=_mm256_hadd_pd(cv,cv);

      res+=rv[1]+rv[2];
    }
    #endif
    #endif

    //add the rest of the vector
    for (;iElement<iElementMax;iElement++) {
      data=ElementDataTable+iElement;
      res+=data->MatrixElementValue*(*data->rhs);
    }

     #if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
     if (irow>=length) exit(__LINE__,__FILE__,"Error: out of bound");
     #endif

     p[irow]=res;
  }
}

}
#endif //_PIC_MATMUL_MPI_MULTITHREAD_
#endif //__CUDA_ARCH__

template <class cCornerNode, int NodeUnknownVariableVectorLength,int MaxStencilLength,
int MaxRhsSupportLength_CornerNodes,int MaxRhsSupportLength_CenterNodes,
int MaxMatrixElementParameterTableLength,int MaxMatrixElementSupportTableLength>
void cLinearSystemCornerNode<cCornerNode, NodeUnknownVariableVectorLength,MaxStencilLength,MaxRhsSupportLength_CornerNodes,MaxRhsSupportLength_CenterNodes,MaxMatrixElementParameterTableLength,MaxMatrixElementSupportTableLength>::Solve(void (*fInitialUnknownValues)(double* x,cCornerNode* CornerNode),
      void (*fUnpackSolution)(double* x,cCornerNode* CornerNode),
      double Tol, //the residual tolerance. The recommended value is 1e-5;
      int nMaxIter, //the max iteration error allowed. The recommended value is 100
      int (*fPackBlockData)(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>** NodeTable,int NodeTableLength,unsigned long int* NodeDataLength,unsigned char* BlockCenterNodeMask,unsigned char* BlockCornerNodeMask,char* SendDataBuffer),
      int (*fUnpackBlockData)(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>** NodeTable,int NodeTableLength,unsigned char* BlockCenterNodeMask,unsigned char* BlockCornerNodeMask,char* RecvDataBuffer)
      ) {
  cMatrixRow* row;
  int cnt;

  //count the total number of the variables, and create the data buffers
  if (SubdomainPartialRHS==NULL) {
    for (SubdomainPartialUnknownsVectorLength=0,row=MatrixRowListFirst;row!=NULL;row=row->next) SubdomainPartialUnknownsVectorLength+=NodeUnknownVariableVectorLength;

    SubdomainPartialRHS=new double [SubdomainPartialUnknownsVectorLength];
    SubdomainPartialUnknownsVector=new double [SubdomainPartialUnknownsVectorLength];
  }

  //calculate the mean value of the Rhs that could be used later for the point located in the 'virtual' blocks
  int i;
  double MeanRhs[NodeUnknownVariableVectorLength];
  int MeanRhsCounter[NodeUnknownVariableVectorLength];

  for (i=0;i<NodeUnknownVariableVectorLength;i++) MeanRhs[i]=0.0,MeanRhsCounter[i]=0;

  for (row=MatrixRowListFirst;row!=NULL;cnt++,row=row->next) if (row->CornerNode!=NULL) {
    MeanRhs[row->iVar]+=row->Rhs;
    MeanRhsCounter[row->iVar]++;
  }

  for (i=0;i<NodeUnknownVariableVectorLength;i++) {
    if (MeanRhsCounter[i]==0) MeanRhs[i]=0.0;
    else MeanRhs[i]/=MeanRhsCounter[i];
  }


  //populate the buffers
  for (cnt=0,row=MatrixRowListFirst;row!=NULL;cnt++,row=row->next) {
    if (row->CornerNode!=NULL) {
      if (cnt%NodeUnknownVariableVectorLength==0) {
        fInitialUnknownValues(SubdomainPartialUnknownsVector+NodeUnknownVariableVectorLength*row->CornerNode->LinearSolverUnknownVectorIndex,row->CornerNode);
      }

      SubdomainPartialRHS[cnt]=row->Rhs;
    }
    else {
      //the point is located in the 'virtual' block
      SubdomainPartialUnknownsVector[row->iVar+NodeUnknownVariableVectorLength*row->ElementDataTable->UnknownVectorIndex]=MeanRhs[row->iVar];
      SubdomainPartialRHS[cnt]=MeanRhs[row->iVar];
    }
  }

  //call the iterative solver
  int nVar=NodeUnknownVariableVectorLength; //variable number
  int nDim = 3; //dimension
  int nI=_BLOCK_CELLS_X_;
  int nJ=_BLOCK_CELLS_Y_;
  int nK=_BLOCK_CELLS_Z_;
  int nBlock = nRealBlocks+nVirtualBlocks;

  MPI_Fint iComm = MPI_Comm_c2f(MPI_GLOBAL_COMMUNICATOR);

  //double Rhs_I(nVar*nVar*nI*nJ*nK*nBlock)// RHS of the equation
  //double Sol_I(nVar*nVar*nI*nJ*nK*nBlock)// vector for solution
  double * Rhs_I=SubdomainPartialRHS;
  double * Sol_I=SubdomainPartialUnknownsVector;
  double PrecondParam=0; // not use preconditioner
//  double ** precond_matrix_II;// pointer to precondition matrix; us  null if no preconditioner
  int lTest=(PIC::ThisThread==0);//1: need test output; 0: no test statement

  linear_solver_wrapper("GMRES", &Tol,&nMaxIter, &nVar, &nDim,&nI, &nJ, &nK, &nBlock, &iComm, Rhs_I,Sol_I, &PrecondParam, NULL, &lTest);

  //unpack the solution
  for (row=MatrixRowListFirst;row!=NULL;row=row->next) if (row->CornerNode!=NULL)  {
    fUnpackSolution(SubdomainPartialUnknownsVector+NodeUnknownVariableVectorLength*row->CornerNode->LinearSolverUnknownVectorIndex,row->CornerNode);
  }

  //execute the data exchange between 'ghost' blocks
  if (fPackBlockData!=NULL) {
    switch (_PIC_BC__PERIODIC_MODE_) {
    case _PIC_BC__PERIODIC_MODE_OFF_:
      PIC::Mesh::mesh->ParallelBlockDataExchange(fPackBlockData,fUnpackBlockData);
      break;

    case _PIC_BC__PERIODIC_MODE_ON_:
      PIC::Parallel::UpdateGhostBlockData(fPackBlockData,fUnpackBlockData);
      break;
    }
  }
  else {
    exit(__LINE__,__FILE__,"Error: fPackBlockData and fUnpackBlockData have to be defined");
  }

}

#endif /* _LINEARSYSTEMCORNERNODE_H_ */
