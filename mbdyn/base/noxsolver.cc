/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2023
 *
 * Pierangelo Masarati  <pierangelo.masarati@polimi.it>
 * Paolo Mantegazza     <paolo.mantegazza@polimi.it>
 *
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 * via La Masa, 34 - 20156 Milano, Italy
 * http://www.aero.polimi.it
 *
 * Changing this copyright notice is forbidden.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation (version 2 of the License).
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/*
  AUTHOR: Reinhard Resch <mbdyn-user@a1.net>
  Copyright (C) 2022(-2023) all rights reserved.

  The copyright of this code is transferred
  to Pierangelo Masarati and Paolo Mantegazza
  for use in the software MBDyn as described
  in the GNU Public License version 2.1

  Tpetra port:
    NOX::Epetra::*              -> NOX::Thyra::*
    AztecOO                     -> Belos (via Stratimikos)
    Epetra_MpiComm/SerialComm   -> Tpetra communicator
    Epetra_Map/Operator/Vector  -> Thyra/Tpetra equivalents

  The linearisation interface (computeF / computeJacobian /
  computePreconditioner) is re-expressed through the NOX::Thyra model
  evaluator thin wrapper ModelEvaluatorWrapper defined below.
*/

#include "mbconfig.h"

#ifdef USE_TRILINOS
#include "solman.h"
#include "solver.h"
#include "noxsolver.h"
#include "output.h"
#include "tpetraspmh.h"
#ifdef USE_MPI
#include "mbcomm.h"
#endif

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wcpp"

#define HAVE_BLAS_SAVE HAVE_BLAS
#define HAVE_BOOL_SAVE HAVE_BOOL
#undef HAVE_BLAS
#undef HAVE_BOOL

/* ---------- Tpetra / Thyra core ---------------------------------------- */
#include "tpetra_types.h"
#include <Thyra_TpetraThyraWrappers.hpp>
#include <Thyra_TpetraLinearOp.hpp>
#include <Thyra_VectorBase.hpp>
#include <Thyra_VectorSpaceBase.hpp>
#include <Thyra_LinearOpWithSolveFactoryBase.hpp>
#include <Thyra_DefaultPreconditioner.hpp>

/* ---------- Stratimikos (Belos via Thyra) ------------------------------- */
#include <Stratimikos_DefaultLinearSolverBuilder.hpp>
#include <Thyra_LinearOpWithSolveFactoryHelpers.hpp>

/* ---------- NOX Thyra interface ----------------------------------------- */
#include <NOX.H>
#include <NOX_Thyra.H>
#include <NOX_Thyra_Group.H>
#include <NOX_Thyra_Vector.H>
#include <NOX_Abstract_PrePostOperator.H>
#include <NOX_Solver_Generic.H>
#include <NOX_Solver_LineSearchBased.H>

/* ---------- Thyra model evaluator --------------------------------------- */
#include <Thyra_ModelEvaluator.hpp>
#include <Thyra_StateFuncModelEvaluatorBase.hpp>

/* ---------- Teuchos ----------------------------------------------------- */
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

#include <BelosOutputManager.hpp>

#undef HAVE_BLAS
#undef HAVE_BOOL
#define HAVE_BLAS HAVE_BLAS_SAVE
#define HAVE_BOOL HAVE_BOOL_SAVE
#undef HAVE_BLAS_SAVE
#undef HAVE_BOOL_SAVE

#pragma GCC diagnostic pop

#ifdef DEBUG_JACOBIAN
#include "sp_gradient_spmh.h"
#endif

/* -------------------------------------------------------------------------
 * NoxSolverParameters
 * ------------------------------------------------------------------------- */
NoxSolverParameters::NoxSolverParameters()
     :CommonNonlinearSolverParam(SOLVER_LINESEARCH_BASED   |
                                 JACOBIAN_NEWTON           |
                                 DIRECTION_NEWTON          |
                                 FORCING_TERM_CONSTANT     |
                                 LINESEARCH_BACKTRACK      |
                                 LINEAR_SOLVER_BLOCK_GMRES |
                                 RECOVERY_STEP_TYPE_CONST,
                                 0,
                                 false),
      dWrmsRelTol(0.),
      dWrmsAbsTol(0.),
      dTolLinSol(1e-10),
      dMinStep(1e-12),
      dRecoveryStep(1.),
      dForcingTermMinTol(1e-6),
      dForcingTermMaxTol(1e-2),
      dForcingTermAlpha(1.5),
      dForcingTermGamma(0.9),
      iMaxIterLinSol(1000),
      iKrylovSubSpaceSize(300),
      iMaxIterLineSearch(200),
      iInnerIterBeforeAssembly(std::numeric_limits<integer>::max())
{
}

namespace {

class NoxNonlinearSolver;

/* =========================================================================
 * Status tests
 * ========================================================================= */
class NoxStatusTest : public NOX::StatusTest::Generic {
public:
     explicit NoxStatusTest(NoxNonlinearSolver& s)
          : oNoxSolver(s), eStatus(NOX::StatusTest::Unevaluated) {}
     virtual ~NoxStatusTest() {}
     NOX::StatusTest::StatusType getStatus() const override { return eStatus; }
     void Reset() { eStatus = NOX::StatusTest::Unevaluated; }
protected:
     NoxNonlinearSolver& oNoxSolver;
     NOX::StatusTest::StatusType eStatus;
};

class NoxResidualTest : public NoxStatusTest {
public:
     explicit NoxResidualTest(NoxNonlinearSolver& s)
          : NoxStatusTest(s), dErrRes(-1.), dErrResDiff(-1.), dTolRes(-2.) {}
     ~NoxResidualTest() {}
     NOX::StatusTest::StatusType
     checkStatus(const NOX::Solver::Generic& problem,
                 NOX::StatusTest::CheckType checkType) override;
     std::ostream& print(std::ostream& stream, int indent) const override;
     void Reset() { NoxStatusTest::Reset(); dErrRes = dErrResDiff = -1.; }
     void SetTolerance(doublereal dTol) { ASSERT(dTol >= 0.); dTolRes = dTol; }
     doublereal dGetTest()     const { return dErrRes;     }
     doublereal dGetTestDiff() const { return dErrResDiff; }
private:
     doublereal dErrRes, dErrResDiff, dTolRes;
};

class NoxSolutionTest : public NoxStatusTest {
public:
     explicit NoxSolutionTest(NoxNonlinearSolver& s)
          : NoxStatusTest(s), dErrSol(-1.), dTolSol(-2.) {}
     ~NoxSolutionTest() {}
     NOX::StatusTest::StatusType
     checkStatus(const NOX::Solver::Generic& problem,
                 NOX::StatusTest::CheckType checkType) override;
     std::ostream& print(std::ostream& stream, int indent) const override;
     void Reset() { NoxStatusTest::Reset(); dErrSol = -1.; }
     void SetTolerance(doublereal dTol) { ASSERT(dTol >= 0.); dTolSol = dTol; }
     doublereal dGetTolerance() const { ASSERT(dTolSol >= 0.); return dTolSol; }
     doublereal dGetTest()      const { return dErrSol; }
private:
     doublereal dErrSol, dTolSol;
};

/* =========================================================================
 * TpetraMatFreeJacOper
 *
 * Matrix-free Jacobian operator for the JACOBIAN_NEWTON_KRYLOV path.
 * Replaces NoxMatrixFreeJacOper (Epetra_Operator) from the original.
 * Implements Thyra::LinearOpBase<SC> so it is handed directly to the
 * NOX::Thyra::Group as the Jacobian and applies J*x on the fly via
 * pNonlinearProblem->Jacobian(&Y, &X).
 * ========================================================================= */
class TpetraMatFreeJacOper : public Thyra::LinearOpBase<TpetraSC>
{
public:
     using SC = TpetraSC;

     TpetraMatFreeJacOper(
          NoxNonlinearSolver& solver,
          const Teuchos::RCP<const Thyra::VectorSpaceBase<SC>>& pSpace_a)
          : oNoxSolver(solver), pSpace(pSpace_a)
#ifdef DEBUG_JACOBIAN
          , pA(nullptr)
#endif
     {}

     ~TpetraMatFreeJacOper()
     {
#ifdef DEBUG_JACOBIAN
          if (pA) { SAFEDELETE(pA); }
#endif
     }

     Teuchos::RCP<const Thyra::VectorSpaceBase<SC>>
     range()  const override { return pSpace; }

     Teuchos::RCP<const Thyra::VectorSpaceBase<SC>>
     domain() const override { return pSpace; }

     bool opSupportedImpl(Thyra::EOpTransp M_trans) const override
     {
          return (M_trans == Thyra::NOTRANS);
     }

protected:
     void applyImpl(
          const Thyra::EOpTransp                          M_trans,
          const Thyra::MultiVectorBase<SC>&               X_in,
          const Teuchos::Ptr<Thyra::MultiVectorBase<SC>>& Y_out,
          const SC                                        alpha,
          const SC                                        beta) const override;

private:
     NoxNonlinearSolver& oNoxSolver;
     Teuchos::RCP<const Thyra::VectorSpaceBase<SC>> pSpace;
#ifdef DEBUG_JACOBIAN
     mutable SpGradientSparseMatrixHandler* pA;
     mutable MyVectorHandler AX;
#endif
};

/* =========================================================================
 * ModelEvaluatorWrapper
 *
 * Bridges MBDyn residual / Jacobian to NOX::Thyra via Thyra::ModelEvaluator.
 *
 * Key design points vs. the Epetra version:
 *
 *  1. create_W_op() returns a Thyra::TpetraLinearOp that shares the same
 *     underlying Tpetra::CrsMatrix that pSolutionManager->pMatHdl() writes
 *     into — assembly in evalModelImpl requires zero data copies.
 *
 *  2. For JACOBIAN_NEWTON_KRYLOV, create_W_op() returns a
 *     TpetraMatFreeJacOper instead.
 *
 *  3. evalModelImpl fully mirrors computeF + computeJacobian from the
 *     original, including:
 *       - residual branch         <- computeF
 *       - W_op branch             <- computeJacobian / Jacobian()
 *       - MBDyn convention that Residual() must precede Jacobian()
 *       - line-search lambda and iteration-counter bookkeeping
 *       - NOX sign convention (residual negated)
 *
 *  4. When MatrInitialize() reallocates the CrsMatrix the cached
 *     pJacobianOp is nulled so create_W_op() re-wraps the new pointer.
 * ========================================================================= */
class ModelEvaluatorWrapper
     : public Thyra::StateFuncModelEvaluatorBase<TpetraSC>
{
public:
     using SC = TpetraSC;

     ModelEvaluatorWrapper(NoxNonlinearSolver& solver,
                           integer iSize,
                           const Teuchos::RCP<const TpetraComm>& pComm_a);

     Teuchos::RCP<const Thyra::VectorSpaceBase<SC>>
     get_x_space() const override { return pSpace; }

     Teuchos::RCP<const Thyra::VectorSpaceBase<SC>>
     get_f_space() const override { return pSpace; }

     Thyra::ModelEvaluatorBase::InArgs<SC>
     getNominalValues() const override { return oNominalValues; }

     Thyra::ModelEvaluatorBase::InArgs<SC>  createInArgs()      const override;
     Thyra::ModelEvaluatorBase::OutArgs<SC> createOutArgsImpl()  const override;

     Teuchos::RCP<Thyra::LinearOpBase<SC>> create_W_op() const override;

     void evalModelImpl(
          const Thyra::ModelEvaluatorBase::InArgs<SC>&  inArgs,
          const Thyra::ModelEvaluatorBase::OutArgs<SC>& outArgs) const override;

     void Rebuild(integer iSize);

     const Teuchos::RCP<const Thyra::VectorSpaceBase<SC>>&
     GetSpace() const { return pSpace; }

private:
     friend class NoxNonlinearSolver; // needs to null pJacobianOp on rebuild

     NoxNonlinearSolver& oNoxSolver;
     Teuchos::RCP<const TpetraComm>                 pComm;
     Teuchos::RCP<const TpetraMap>                  pMap;
     Teuchos::RCP<const Thyra::VectorSpaceBase<SC>> pSpace;
     Thyra::ModelEvaluatorBase::InArgs<SC>          oNominalValues;

     // Lazily initialised; mutable because create_W_op() is const.
     // Nulled by Rebuild() and by Jacobian() after ErrRebuildMatrix.
     mutable Teuchos::RCP<Thyra::LinearOpBase<SC>>  pJacobianOp;
};

/* =========================================================================
 * NoxNonlinearSolver
 * ========================================================================= */
class NoxNonlinearSolver : public NonlinearSolver,
                           private NoxSolverParameters,
                           private NOX::Abstract::PrePostOperator
{
public:
     friend class ModelEvaluatorWrapper;
     friend class TpetraMatFreeJacOper;
     friend class NoxResidualTest;
     friend class NoxSolutionTest;

     NoxNonlinearSolver(const NonlinearSolverTestOptions& oSolverOpt,
                        const NoxSolverParameters& oParam);
     ~NoxNonlinearSolver();

     void Solve(const NonlinearProblem* pNLP,
                Solver* pS,
                const integer iMaxIter,
                const doublereal& Tol,
                integer& iIterCnt,
                doublereal& dErr,
                const doublereal& SolTol,
                doublereal& dSolErr) override;

     bool NoxMakeSolTest(const VectorHandler& XPrev,
                         const VectorHandler& XCurr,
                         const doublereal& dTol,
                         doublereal& dTest);

     bool NoxMakeResTest(const VectorHandler& oResVec,
                         const doublereal& dTol,
                         doublereal& dTest,
                         doublereal& dTestDiff);

private:
     struct CPUTimeGuard {
          explicit CPUTimeGuard(const NoxNonlinearSolver& s, CPUTimeType t)
               : oWatch(s, t) { oWatch.Tic(); }
          ~CPUTimeGuard() { oWatch.Toc(); }
          CPUStopWatch oWatch;
     };

     void Attach(Solver* pS, const NonlinearProblem* pNLP);
     void Residual(const VectorHandler* pSol, VectorHandler* pRes);
     void Jacobian();
     void BuildSolver(integer iMaxIter);
     void OutputIteration(integer iIterCnt, bool bJacobian) const;

     void runPreIterate(const NOX::Solver::Generic&)            override {}
     void runPostIterate(const NOX::Solver::Generic&)           override {}
     void runPreSolve(const NOX::Solver::Generic&)              override {}
     void runPostSolve(const NOX::Solver::Generic&)             override {}
     void runPreSolutionUpdate(const NOX::Abstract::Vector&,
                               const NOX::Solver::Generic&)     override {}
     void runPreLineSearch(const NOX::Solver::Generic& solver)  override;
     void runPostLineSearch(const NOX::Solver::Generic& solver) override;

     void ResetPrecondReuse() const {
          iPrecInnerIterCnt = iPrecInnerIterCntTot = 0;
     }
     void ForcePrecondRebuild() const {
          iPrecInnerIterCnt = iPrecInnerIterCntTot =
               std::numeric_limits<integer>::max();
     }

     Teuchos::RCP<const TpetraComm>      pComm;
     Teuchos::RCP<ModelEvaluatorWrapper> pModelEval;
     Teuchos::RCP<NOX::Solver::Generic>  pNonlinearSolver;
     Teuchos::RCP<NOX::Thyra::Vector>    pSolutionView;
     Teuchos::ParameterList              oSolverParam;

     const NonlinearProblem* pNonlinearProblem;
     Solver*                 pSolver;
     SolutionManager*        pSolutionManager;

     mutable MyVectorHandler DeltaX, XPrev, TmpRes;

     NoxResidualTest oResTest;
     NoxSolutionTest oSolTest;

     mutable bool bUseTranspose;
     bool         bUpdateJacobian;
     bool         bInDerivativeSolver;
     bool         bInLineSearch;

     mutable integer iPrecInnerIterCnt;
     mutable integer iPrecInnerIterCntTot;
     mutable integer iInnerIterCntTot;
};

/* =========================================================================
 * TpetraMatFreeJacOper::applyImpl
 * ========================================================================= */
void TpetraMatFreeJacOper::applyImpl(
     const Thyra::EOpTransp                          M_trans,
     const Thyra::MultiVectorBase<SC>&               X_in,
     const Teuchos::Ptr<Thyra::MultiVectorBase<SC>>& Y_out,
     const SC                                        alpha,
     const SC                                        beta) const
{
     if (M_trans != Thyra::NOTRANS) {
          throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
     }
     // NOX always calls with alpha=1, beta=0
     ASSERT(alpha == SC(1.) && beta == SC(0.));

     NoxNonlinearSolver::CPUTimeGuard oCPUTimeJac(
          oNoxSolver, NoxNonlinearSolver::CPU_JACOBIAN);

     ASSERT(oNoxSolver.Size > 0);

     auto xTpetra =
          Thyra::TpetraOperatorVectorExtraction<SC, TpetraLO, TpetraGO, TpetraNode>
               ::getConstTpetraMultiVector(Teuchos::rcpFromRef(X_in));
     auto yTpetra =
          Thyra::TpetraOperatorVectorExtraction<SC, TpetraLO, TpetraGO, TpetraNode>
               ::getTpetraMultiVector(Teuchos::rcpFromRef(*Y_out));

     ASSERT(static_cast<integer>(xTpetra->getLocalLength()) == oNoxSolver.Size);
     ASSERT(static_cast<integer>(yTpetra->getLocalLength()) == oNoxSolver.Size);

     const MyVectorHandler X(
          oNoxSolver.Size,
          const_cast<doublereal*>(
               xTpetra->getLocalViewHost(Tpetra::Access::ReadOnly).data()));
     MyVectorHandler Y(
          oNoxSolver.Size,
          yTpetra->getLocalViewHost(Tpetra::Access::ReadWrite).data());

     // Matrix-free product: Y = J * X
     oNoxSolver.pNonlinearProblem->Jacobian(&Y, &X);

#ifdef DEBUG_JACOBIAN
     if (!pA) {
          SAFENEWWITHCONSTRUCTOR(pA,
                                  SpGradientSparseMatrixHandler,
                                  SpGradientSparseMatrixHandler(
                                       oNoxSolver.Size, oNoxSolver.Size));
          oNoxSolver.pNonlinearProblem->Jacobian(pA);
          pA->PacMat();
          AX.Resize(oNoxSolver.Size);
     }
     pA->MatVecMul(AX, X);
     const doublereal dTolRel = sqrt(std::numeric_limits<doublereal>::epsilon());
     const doublereal dTolAbs = sqrt(std::numeric_limits<doublereal>::epsilon());
     const doublereal dNormAX = AX.Norm();
     for (integer i = 1; i <= X.iGetSize(); ++i) {
          if (std::fabs(AX(i) - Y(i)) > dTolRel + dTolAbs * dNormAX) {
               DEBUGCERR("MatFree Jacobian check failed: AX("
                         << i << ")=" << AX(i)
                         << " Y(" << i << ")=" << Y(i) << "\n");
               ASSERT(0);
          }
     }
#endif
}

/* =========================================================================
 * ModelEvaluatorWrapper  implementation
 * ========================================================================= */
ModelEvaluatorWrapper::ModelEvaluatorWrapper(
     NoxNonlinearSolver& solver,
     integer iSize,
     const Teuchos::RCP<const TpetraComm>& pComm_a)
     : oNoxSolver(solver), pComm(pComm_a)
{
     Rebuild(iSize);
}

void ModelEvaluatorWrapper::Rebuild(integer iSize)
{
     pMap = Teuchos::rcp(new TpetraMap(
          static_cast<TpetraGO>(iSize), TpetraGO(0), pComm));
     pSpace =
          Thyra::createVectorSpace<TpetraSC, TpetraLO, TpetraGO, TpetraNode>(pMap);

     pJacobianOp = Teuchos::null; // invalidate on resize

     oNominalValues = this->createInArgs();
     auto x0 = Thyra::createMember(pSpace);
     Thyra::assign(x0.ptr(), TpetraSC(0.));
     oNominalValues.set_x(x0);
}

Thyra::ModelEvaluatorBase::InArgs<TpetraSC>
ModelEvaluatorWrapper::createInArgs() const
{
     Thyra::ModelEvaluatorBase::InArgsSetup<SC> inArgs;
     inArgs.setModelEvalDescription(this->description());
     inArgs.setSupports(Thyra::ModelEvaluatorBase::IN_ARG_x);
     return inArgs;
}

Thyra::ModelEvaluatorBase::OutArgs<TpetraSC>
ModelEvaluatorWrapper::createOutArgsImpl() const
{
     Thyra::ModelEvaluatorBase::OutArgsSetup<SC> outArgs;
     outArgs.setModelEvalDescription(this->description());
     outArgs.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_f);
     outArgs.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_W_op);
     return outArgs;
}

Teuchos::RCP<Thyra::LinearOpBase<TpetraSC>>
ModelEvaluatorWrapper::create_W_op() const
{
     // Lazily construct and cache.
     // pSolutionManager must already be valid (Attach() precedes BuildSolver()).
     if (pJacobianOp.is_null()) {
          if (oNoxSolver.uFlags & NoxSolverParameters::JACOBIAN_NEWTON_KRYLOV) {
               // Matrix-free: recomputes J*v on every apply call.
               pJacobianOp = Teuchos::rcp(
                    new TpetraMatFreeJacOper(oNoxSolver, pSpace));
          } else {
               // Explicit matrix: wrap the Tpetra::CrsMatrix that the solution
               // manager owns — shared pointer means zero-copy assembly.
               ASSERT(oNoxSolver.pSolutionManager != nullptr);
               MatrixHandler* pMH = oNoxSolver.pSolutionManager->pMatHdl();
               auto* pTpetraMH = dynamic_cast<TpetraSparseMatrixHandler*>(pMH);
               ASSERT(pTpetraMH != nullptr);
               pJacobianOp =
                    Thyra::tpetraLinearOp<TpetraSC, TpetraLO, TpetraGO, TpetraNode>(
                         pSpace, pSpace,
                         pTpetraMH->pGetTpetraCrsMatrix()); // adapt to actual API
          }
     }
     return pJacobianOp;
}

void ModelEvaluatorWrapper::evalModelImpl(
     const Thyra::ModelEvaluatorBase::InArgs<SC>&  inArgs,
     const Thyra::ModelEvaluatorBase::OutArgs<SC>& outArgs) const
{
     /* --- Unpack x (zero-copy) ------------------------------------------ */
     auto xThyra = inArgs.get_x();
     ASSERT(!xThyra.is_null());

     auto xTpetra =
          Thyra::TpetraOperatorVectorExtraction<SC, TpetraLO, TpetraGO, TpetraNode>
               ::getConstTpetraVector(xThyra);

     const MyVectorHandler oSol(
          oNoxSolver.Size,
          const_cast<doublereal*>(
               xTpetra->getLocalViewHost(Tpetra::Access::ReadOnly).data()));

     /* --- Residual branch  (mirrors computeF) --------------------------- */
     if (!outArgs.get_f().is_null()) {
          auto fThyra  = outArgs.get_f();
          auto fTpetra =
               Thyra::TpetraOperatorVectorExtraction<SC, TpetraLO, TpetraGO, TpetraNode>
                    ::getTpetraVector(fThyra);

          MyVectorHandler oRes(
               oNoxSolver.Size,
               fTpetra->getLocalViewHost(Tpetra::Access::ReadWrite).data());

          oNoxSolver.bUpdateJacobian = true;

          if (oNoxSolver.bInLineSearch && oNoxSolver.pNonlinearSolver) {
               const auto& oLS =
                    dynamic_cast<const NOX::Solver::LineSearchBased&>(
                         *oNoxSolver.pNonlinearSolver);
               oNoxSolver.SetNonlinearSolverHint(
                    NonlinearSolver::LINESEARCH_LAMBDA_CURR,
                    oLS.getStepSize());
               DEBUGCERR("line search iteration "
                    << oNoxSolver.GetNonlinearSolverHint(
                         NonlinearSolver::LINESEARCH_ITERATION_CURR)
                    << ": lambda=" << oLS.getStepSize() << "\n");
          }

          oNoxSolver.Residual(&oSol, &oRes);
          // std::cerr << "chiamato Residual" << std::endl;
          //           std::cerr << "==============" << std::endl;
          //           oNoxSolver.pSolver->PrintSolution(
          //                oNoxSolver.DeltaX,
          //                oNoxSolver.pNonlinearSolver->getNumIterations());
          //           std::cerr << "--------------" << std::endl;
          //           std::cerr << oSol << std::endl;
          //           std::cerr << ";;;;;;;;;;;;;;" << std::endl;
          //           oNoxSolver.pSolver->PrintResidual(
          //                oRes,
          //                oNoxSolver.pNonlinearSolver->getNumIterations());
          //           std::cerr << "**************" << std::endl;
          //           std::cerr << oRes << std::endl;
          //           std::cerr << "++++++++++++++" << std::endl;

          if (oNoxSolver.pSolver && oNoxSolver.pNonlinearSolver) {
               if (oNoxSolver.outputSol()) {
                    oNoxSolver.pSolver->PrintSolution(
                         oNoxSolver.DeltaX,
                         oNoxSolver.pNonlinearSolver->getNumIterations());
               }
               if (oNoxSolver.outputRes()) {
                    oNoxSolver.pSolver->PrintResidual(
                         oRes,
                         oNoxSolver.pNonlinearSolver->getNumIterations());
               }
          }

          oRes *= -1.; // NOX sign convention

          if (oNoxSolver.bInLineSearch) {
               const integer iIterCurr =
                    oNoxSolver.GetNonlinearSolverHint(
                         NonlinearSolver::LINESEARCH_ITERATION_CURR);
               oNoxSolver.SetNonlinearSolverHint(
                    NonlinearSolver::LINESEARCH_ITERATION_CURR,
                    iIterCurr + 1);
          }
     }

     /* --- Jacobian branch  (mirrors computeJacobian) -------------------- */
     // outArgs.get_W_op() returns the same Thyra::TpetraLinearOp that
     // create_W_op() gave the NOX::Thyra::Group.  Because it wraps the
     // solution manager's CrsMatrix by shared_ptr, assembling below into
     // pSolutionManager->pMatHdl() automatically updates the operator NOX
     // uses for the linear solve — no copy required.
     // For the matrix-free path no matrix assembly is performed at all.
     if (!outArgs.get_W_op().is_null()) {
          if (oNoxSolver.bUpdateJacobian) {
               // MBDyn convention: Residual() must precede Jacobian().
               // If the residual branch above already ran, the solution is
               // current; otherwise call explicitly (preconditioner refresh).
               if (outArgs.get_f().is_null()) {
                    oNoxSolver.Residual(&oSol, &oNoxSolver.TmpRes);
               }
               if (oNoxSolver.pSolutionManager) {
                    oNoxSolver.pSolutionManager->MatrReset();
               }
               oNoxSolver.Jacobian(); // writes into pSolutionManager->pMatHdl()
               oNoxSolver.bUpdateJacobian = false;
          }
     }
}

/* =========================================================================
 * NoxResidualTest  implementation
 * ========================================================================= */
NOX::StatusTest::StatusType
NoxResidualTest::checkStatus(const NOX::Solver::Generic& problem,
                              NOX::StatusTest::CheckType  checkType)
{
     if (checkType == NOX::StatusTest::None) {
          eStatus = NOX::StatusTest::Unevaluated;
          return eStatus;
     }
     const NOX::Abstract::Group& grp = problem.getSolutionGroup();
     if (!grp.isF()) {
          eStatus = NOX::StatusTest::Unevaluated;
          return eStatus;
     }

     const auto& FT =
          dynamic_cast<const NOX::Thyra::Vector&>(grp.getF());
     auto fTpetra =
          Thyra::TpetraOperatorVectorExtraction<TpetraSC, TpetraLO, TpetraGO, TpetraNode>
               ::getConstTpetraVector(
                    Teuchos::rcpFromRef(FT.getThyraVector()));

     const MyVectorHandler oResVec(
          static_cast<integer>(fTpetra->getLocalLength()),
          const_cast<doublereal*>(
               fTpetra->getLocalViewHost(Tpetra::Access::ReadOnly).data()));

     const doublereal dFirstResFact =
          problem.getNumIterations() == 0 ? 1e-2 : 1.;

     eStatus = oNoxSolver.NoxMakeResTest(
          oResVec, dFirstResFact * dTolRes, dErrRes, dErrResDiff)
          ? NOX::StatusTest::Converged
          : NOX::StatusTest::Unconverged;
     return eStatus;
}

std::ostream& NoxResidualTest::print(std::ostream& stream, int indent) const
{
     for (int j = 0; j < indent; ++j) stream << ' ';
     stream << eStatus
            << "F-Norm = " << NOX::Utils::sciformat(dErrRes, 3)
            << " < "       << NOX::Utils::sciformat(dTolRes, 3) << "\n";
     return stream;
}

/* =========================================================================
 * NoxSolutionTest  implementation
 * ========================================================================= */
NOX::StatusTest::StatusType
NoxSolutionTest::checkStatus(const NOX::Solver::Generic& problem,
                              NOX::StatusTest::CheckType  checkType)
{
     if (checkType == NOX::StatusTest::None) {
          eStatus = NOX::StatusTest::Unevaluated;
          return eStatus;
     }
     if (problem.getNumIterations() == 0) {
          eStatus = NOX::StatusTest::Unconverged;
          return eStatus;
     }
     if (!problem.getSolutionGroup().isF()) {
          eStatus = NOX::StatusTest::Unconverged;
          return eStatus;
     }

     const integer n = oNoxSolver.Size;
     auto extractRaw = [n](const NOX::Abstract::Vector& VA) {
          const auto& VT = dynamic_cast<const NOX::Thyra::Vector&>(VA);
          auto tp =
               Thyra::TpetraOperatorVectorExtraction<
                    TpetraSC, TpetraLO, TpetraGO, TpetraNode>
                         ::getConstTpetraVector(
                              Teuchos::rcpFromRef(VT.getThyraVector()));
          return MyVectorHandler(
               n,
               const_cast<doublereal*>(
                    tp->getLocalViewHost(Tpetra::Access::ReadOnly).data()));
     };

     MyVectorHandler XPrevVec =
          extractRaw(problem.getPreviousSolutionGroup().getX());
     MyVectorHandler XCurrVec =
          extractRaw(problem.getSolutionGroup().getX());

     eStatus = oNoxSolver.NoxMakeSolTest(XPrevVec, XCurrVec, dTolSol, dErrSol)
          ? NOX::StatusTest::Converged
          : NOX::StatusTest::Unconverged;
     return eStatus;
}

std::ostream& NoxSolutionTest::print(std::ostream& stream, int indent) const
{
     for (int j = 0; j < indent; ++j) stream << ' ';
     stream << eStatus
            << "X-Norm = " << NOX::Utils::sciformat(dErrSol, 3)
            << " < "       << NOX::Utils::sciformat(dTolSol, 3) << "\n";
     return stream;
}

/* =========================================================================
 * NoxNonlinearSolver  implementation
 * ========================================================================= */
NoxNonlinearSolver::NoxNonlinearSolver(
     const NonlinearSolverTestOptions& oSolverOpt,
     const NoxSolverParameters& oParam)
     : NonlinearSolver(oSolverOpt),
       NoxSolverParameters(oParam),
       pNonlinearProblem(nullptr),
       pSolver(nullptr),
       pSolutionManager(nullptr),
       oResTest(*this),
       oSolTest(*this),
       bUseTranspose(false),
       bUpdateJacobian(true),
       bInDerivativeSolver(true),
       bInLineSearch(false),
       iInnerIterCntTot(0)
{
#ifdef USE_MPI
     pComm = tpetraMpiComm(MBDynComm);
#else
     pComm = tpetraSerialComm();
#endif
     ForcePrecondRebuild();
}

NoxNonlinearSolver::~NoxNonlinearSolver()
{
     silent_cerr("total inner iterations: " << iInnerIterCntTot << "\n");
}

void
NoxNonlinearSolver::Solve(const NonlinearProblem* pNLP,
                           Solver* pS,
                           const integer iMaxIter,
                           const doublereal& dTolRes,
                           integer& iIterCnt,
                           doublereal& dResErr,
                           const doublereal& dTolSol,
                           doublereal& dSolErr)
{
     DEBUGCERR("Solve()\n");

     bInLineSearch = false;
     SetNonlinearSolverHint(LINESEARCH_ITERATION_CURR, 0);
     SetNonlinearSolverHint(LINESEARCH_LAMBDA_CURR, 1.);

     oResTest.SetTolerance(dTolRes);
     oSolTest.SetTolerance(dTolSol);

     Attach(pS, pNLP);

     if (!pNonlinearSolver.get()) {
          BuildSolver(iMaxIter);
     }

     pNonlinearSolver->reset(*pSolutionView);

     if (!bKeepJacAcrossSteps || bInDerivativeSolver) {
          ForcePrecondRebuild();
     }

     iIterCnt = 0;

     for (;;) {
          oResTest.Reset();
          oSolTest.Reset();

          const integer TotJacPrev = TotJac;

          DEBUGCERR("Nonlinear solver step(" << iIterCnt << ")\n");

          NOX::StatusTest::StatusType solvStatus = pNonlinearSolver->step();

          if (oResTest.getStatus() == NOX::StatusTest::Unevaluated) {
               oResTest.checkStatus(*pNonlinearSolver,
                                    NOX::StatusTest::CheckType::Complete);
          }

          ++iIterCnt;

          OutputIteration(iIterCnt, TotJac > TotJacPrev);
          pSolver->CheckTimeStepLimit(oResTest.dGetTest(),
                                       oResTest.dGetTestDiff());

          dResErr = oResTest.dGetTest();
          dSolErr = (oSolTest.getStatus() != NOX::StatusTest::Unevaluated)
               ? oSolTest.dGetTest() : -1.;

          if (solvStatus == NOX::StatusTest::Converged) break;

          if (iIterCnt > iMaxIter || solvStatus == NOX::StatusTest::Failed) {
               throw NoConvergence(MBDYN_EXCEPT_ARGS);
          }

          if (mbdyn_stop_at_end_of_iteration()) {
               throw ErrInterrupted(MBDYN_EXCEPT_ARGS);
          }
     }
}

bool NoxNonlinearSolver::NoxMakeSolTest(const VectorHandler& XPrev_a,
                                         const VectorHandler& XCurr,
                                         const doublereal& dTol,
                                         doublereal& dTest)
{
     DeltaX.ScalarAddMul(XCurr, XPrev_a, -1.);
     return NonlinearSolver::MakeSolTest(pSolver, DeltaX, dTol, dTest);
}

bool NoxNonlinearSolver::NoxMakeResTest(const VectorHandler& oResVec,
                                         const doublereal& dTol,
                                         doublereal& dTest,
                                         doublereal& dTestDiff)
{
     return NonlinearSolver::MakeResTest(pSolver, pNonlinearProblem,
                                          oResVec, dTol, dTest, dTestDiff);
}

void NoxNonlinearSolver::BuildSolver(const integer iMaxIter_a)
{
     // Attach() must have been called first so pSolutionManager is valid
     // when create_W_op() resolves the Tpetra::CrsMatrix pointer.
     if (!pModelEval) {
          pModelEval = Teuchos::rcp(
               new ModelEvaluatorWrapper(*this, Size, pComm));
     } else {
          pModelEval->Rebuild(Size);
     }

     auto x0 = Thyra::createMember(pModelEval->GetSpace());
     Thyra::assign(x0.ptr(), TpetraSC(0.));
     pSolutionView = Teuchos::rcp(new NOX::Thyra::Vector(x0));

     /* ---- Stratimikos / Belos ------------------------------------------ */
     Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;
     Teuchos::RCP<Teuchos::ParameterList> pLSParams =
          Teuchos::rcp(new Teuchos::ParameterList());

     pLSParams->set("Linear Solver Type", "Belos");

     std::string sBelosSolverType = "Block GMRES";
     if      (uFlags & LINEAR_SOLVER_PSEUDO_BLOCK_GMRES)         sBelosSolverType = "Pseudo Block GMRES";
     else if (uFlags & LINEAR_SOLVER_BLOCK_CG)                   sBelosSolverType = "Block CG";
     else if (uFlags & LINEAR_SOLVER_PSEUDO_BLOCK_CG)            sBelosSolverType = "Pseudo Block CG";
     else if (uFlags & LINEAR_SOLVER_BLOCK_STOCHASTIC_CG)        sBelosSolverType = "Block Stochastic CG";
     else if (uFlags & LINEAR_SOLVER_GCRODR)                     sBelosSolverType = "GCRODR";
     else if (uFlags & LINEAR_SOLVER_RCG)                        sBelosSolverType = "RCG";
     else if (uFlags & LINEAR_SOLVER_MINRES)                     sBelosSolverType = "MINRES";
     else if (uFlags & LINEAR_SOLVER_TFQMR)                      sBelosSolverType = "TFQMR";
     else if (uFlags & LINEAR_SOLVER_BICGSTAB)                   sBelosSolverType = "BiCGStab";
     else if (uFlags & LINEAR_SOLVER_FIXED_POINT)                sBelosSolverType = "Fixed Point";
     else if (uFlags & LINEAR_SOLVER_TPETRA_GMRES)               sBelosSolverType = "TPETRA GMRES";
     else if (uFlags & LINEAR_SOLVER_TPETRA_GMRES_PIPELINE)      sBelosSolverType = "TPETRA GMRES PIPELINE";
     else if (uFlags & LINEAR_SOLVER_TPETRA_GMRES_SINGLE_REDUCE) sBelosSolverType = "TPETRA GMRES SINGLE REDUCE";
     else if (uFlags & LINEAR_SOLVER_TPETRA_GMRES_SSTEP)         sBelosSolverType = "TPETRA GMRES S-STEP";

     const integer iKrylovRestart =
          std::min(Size, std::min(iMaxIterLinSol, iKrylovSubSpaceSize));

     auto& belosParams =
          pLSParams->sublist("Linear Solver Types").sublist("Belos");
     belosParams.set("Solver Type", sBelosSolverType);
     auto& bSolverParams =
          belosParams.sublist("Solver Types").sublist(sBelosSolverType);
     bSolverParams.set("Maximum Iterations",    iMaxIterLinSol);
     bSolverParams.set("Convergence Tolerance", dTolLinSol);
     bSolverParams.set("Num Blocks",            iKrylovRestart);
     bSolverParams.set("Output Frequency",
                       (uFlags & PRINT_CONVERGENCE_INFO) ? 1 : 0);
     bSolverParams.set("Output Style", 1);
     bSolverParams.set("Verbosity",
                       (uFlags & PRINT_CONVERGENCE_INFO)
                            ? (Belos::Errors | Belos::Warnings |
                               Belos::IterationDetails |
                               Belos::StatusTestDetails)
                            : Belos::Errors);

     if (!(uFlags & JACOBIAN_NEWTON_KRYLOV)) {
          pLSParams->set("Preconditioner Type", "Ifpack2");
          pLSParams->sublist("Preconditioner Types")
                    .sublist("Ifpack2")
                    .set("Prec Type", "RILUK")
                    .sublist("Ifpack2 Settings")
                    .set("fact: iluk level-of-fill", 1);
     } else {
          pLSParams->set("Preconditioner Type", "None");
     }

     linearSolverBuilder.setParameterList(pLSParams);
     auto pLOWSFactory = linearSolverBuilder.createLinearSolveStrategy("");

     /* ---- NOX parameters ----------------------------------------------- */
     auto& oPrintParam    = oSolverParam.sublist("Printing");
     auto& oDirectionParam = oSolverParam.sublist("Direction");
     auto& oNewtonParam    = oDirectionParam.sublist("Newton");

     oNewtonParam.set("Forcing Term Minimum Tolerance", dForcingTermMinTol);
     oNewtonParam.set("Forcing Term Maximum Tolerance", dForcingTermMaxTol);
     oNewtonParam.set("Forcing Term Alpha",              dForcingTermAlpha);
     oNewtonParam.set("Forcing Term Gamma",              dForcingTermGamma);

     int iSolverOutput = 0;
     if (outputIters()) {
          if (uFlags & VERBOSE_MODE)           iSolverOutput |= NOX::Utils::Warning;
          if (uFlags & PRINT_CONVERGENCE_INFO) iSolverOutput |=
               NOX::Utils::OuterIteration | NOX::Utils::InnerIteration |
               NOX::Utils::LinearSolverDetails | NOX::Utils::StepperIteration |
               NOX::Utils::StepperDetails | NOX::Utils::Parameters |
               NOX::Utils::Details | NOX::Utils::OuterIterationStatusTest |
               NOX::Utils::TestDetails;
     }
     oPrintParam.set("Output Information", iSolverOutput);

     NOX::Abstract::PrePostOperator& oPrePost = *this;
     oSolverParam.sublist("Solver Options")
          .set("User Defined Pre/Post Operator", Teuchos::rcpFromRef(oPrePost));
     if (bInDerivativeSolver) {
          oSolverParam.sublist("Solver Options")
               .set("Status Test Check Type", "Complete");
     }

     static constexpr char szNLSolver[] = "Nonlinear Solver";

     if (uFlags & SOLVER_LINESEARCH_BASED) {
          oSolverParam.set(szNLSolver, "Line Search Based");
          auto& oLSP = oSolverParam.sublist("Line Search");
          std::string strLSMethod;
          if      (uFlags & LINESEARCH_BACKTRACK)    strLSMethod = "Backtrack";
          else if (uFlags & LINESEARCH_POLYNOMIAL)   strLSMethod = "Polynomial";
          else if (uFlags & LINESEARCH_MORE_THUENTE) strLSMethod = "More'-Thuente";
          else                                       strLSMethod = "Full Step";
          oLSP.set("Method", strLSMethod);
          if      (uFlags & SUFFICIENT_DEC_COND_ARMIJO_GOLDSTEIN)
               oLSP.sublist(strLSMethod).set("Sufficient Decrease Condition","Armijo-Goldstein");
          else if (uFlags & SUFFICIENT_DEC_COND_ARED_PRED)
               oLSP.sublist(strLSMethod).set("Sufficient Decrease Condition","Ared/Pred");
          auto& oLSM = oLSP.sublist(strLSMethod);
          oLSM.set("Max Iters",     iMaxIterLineSearch);
          oLSM.set("Minimum Step",  dMinStep);
          oLSM.set("Recovery Step", dRecoveryStep);
          if      (uFlags & RECOVERY_STEP_TYPE_CONST)
               oLSM.set("Recovery Step Type","Constant");
          else if (uFlags & RECOVERY_STEP_TYPE_LAST_STEP)
               oLSM.set("Recovery Step Type","Last Computed Step");
     } else if (uFlags & SOLVER_TRUST_REGION_BASED) {
          oSolverParam.set(szNLSolver, "Trust Region Based");
     } else if (uFlags & SOLVER_INEXACT_TRUST_REGION_BASED) {
          oSolverParam.set(szNLSolver, "Inexact Trust Region Based");
     } else if (uFlags & SOLVER_TENSOR_BASED) {
          oSolverParam.set(szNLSolver, "Tensor Based");
          oSolverParam.sublist("Line Search").set("Method","Curvilinear")
               .sublist("Curvilinear").set("Minimum Step", dMinStep);
     }

     if (uFlags & DIRECTION_NEWTON) {
          oDirectionParam.set("Method", "Newton");
          if      (uFlags & FORCING_TERM_CONSTANT) oNewtonParam.set("Forcing Term Method","Constant");
          else if (uFlags & FORCING_TERM_TYPE1)    oNewtonParam.set("Forcing Term Method","Type 1");
          else if (uFlags & FORCING_TERM_TYPE2)    oNewtonParam.set("Forcing Term Method","Type 2");
     } else if (uFlags & DIRECTION_STEEPEST_DESCENT) {
          oDirectionParam.set("Method", "Steepest Descent");
     } else if (uFlags & DIRECTION_NONLINEAR_CG) {
          oDirectionParam.set("Method", "NonlinearCG");
     } else if (uFlags & DIRECTION_BROYDEN) {
          oDirectionParam.set("Method", "Broyden");
          auto& oBroyden = oDirectionParam.sublist("Broyden");
          oBroyden.set("Restart Frequency", iIterationsBeforeAssembly);
          if      (uFlags & FORCING_TERM_CONSTANT) oBroyden.set("Forcing Term Method","Constant");
          else if (uFlags & FORCING_TERM_TYPE1)    oBroyden.set("Forcing Term Method","Type 1");
          else if (uFlags & FORCING_TERM_TYPE2)    oBroyden.set("Forcing Term Method","Type 2");
     }

     /* ---- NOX::Thyra::Group -------------------------------------------- */
     Teuchos::RCP<const Thyra::ModelEvaluator<TpetraSC>> pModelEvalConst = pModelEval;
     Teuchos::RCP<Thyra::LinearOpBase<TpetraSC>>         pJacobian = pModelEval->create_W_op();
     Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<TpetraSC>>
          pLOWSFactoryConst = pLOWSFactory;
     Teuchos::RCP<Thyra::PreconditionerBase<TpetraSC>>        pNullPrecOp;
     Teuchos::RCP<Thyra::PreconditionerFactoryBase<TpetraSC>> pNullPrecFactory;

     auto grpPtr = Teuchos::rcp(new NOX::Thyra::Group(
          *pSolutionView, pModelEvalConst, pJacobian,
          pLOWSFactoryConst, pNullPrecOp, pNullPrecFactory));

     /* ---- Status tests -------------------------------------------------- */
     auto converged =
          Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND));
     converged->addStatusTest(Teuchos::rcpFromRef(oResTest));
     if (oSolTest.dGetTolerance() > 0.) {
          converged->addStatusTest(Teuchos::rcpFromRef(oSolTest));
     }
     if (dWrmsRelTol > 0. && dWrmsAbsTol > 0. && !bInDerivativeSolver) {
          converged->addStatusTest(
               Teuchos::rcp(new NOX::StatusTest::NormWRMS(dWrmsRelTol, dWrmsAbsTol)));
     }

     auto pCombCriteria =
          Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
     pCombCriteria->addStatusTest(Teuchos::rcp(new NOX::StatusTest::FiniteValue));
     pCombCriteria->addStatusTest(converged);
     pCombCriteria->addStatusTest(
          Teuchos::rcp(new NOX::StatusTest::MaxIters(iMaxIter_a)));

     ForcePrecondRebuild();
     pNonlinearSolver = NOX::Solver::buildSolver(
          grpPtr, pCombCriteria, Teuchos::rcpFromRef(oSolverParam));
}

void NoxNonlinearSolver::OutputIteration(integer iIterCnt, bool bJacobian) const
{
#ifdef USE_MPI
     if (!bParallel || MBDynComm.Get_rank() == 0)
#endif
     {
          if (!(outputIters() || outputSolverConditionNumber())) return;
          silent_cout("\tIteration(" << iIterCnt << ") " << oResTest.dGetTest());
          if (bJacobian) {
               silent_cout(" J");
               if (outputSolverConditionNumber()) {
                    silent_cout(" cond=");
                    doublereal dCond;
                    if (pSolutionManager->bGetConditionNumber(dCond)) {
                         silent_cout(dCond);
                         if (outputSolverConditionStat()) {
                              AddCond(dCond);
                              silent_cout(" " << dGetCondMin()
                                         << " " << dGetCondMax()
                                         << " " << dGetCondAvg());
                         }
                    } else { silent_cout("NA"); }
               }
          }
          if (outputCPUTime()) {
               typedef std::chrono::duration<float,std::ratio<1,1>> FloatSec;
               auto flags = std::cout.flags(); auto prec = std::cout.precision();
               std::cout.setf(std::ios::scientific); std::cout.precision(2);
               silent_cout(" CPU:"
                    << FloatSec(dGetTimeCPU(CPU_RESIDUAL)).count()   << '+'
                    << FloatSec(dGetTimeCPU(CPU_JACOBIAN)).count()   << '+'
                    << FloatSec(dGetTimeCPU(CPU_LINEAR_SOLVER)).count());
               std::cout.flags(flags); std::cout.precision(prec);
          }
          if (oSolTest.getStatus() != NOX::StatusTest::Unevaluated) {
               silent_cout("\n\t\tSolErr " << oSolTest.dGetTest());
          }
          silent_cout('\n');
     }
}

void NoxNonlinearSolver::Attach(Solver* pS, const NonlinearProblem* pNLP)
{
     pSolver = pS;
     pSolutionManager = pSolver->pGetSolutionManager();
     if (pNLP != pNonlinearProblem) {
          DEBUGCERR("Resetting nonlinear solver\n");
          pNonlinearSolver.reset();
          bUpdateJacobian = true;
          ResetCond();
          if (pNonlinearProblem) bInDerivativeSolver = false;
     }
     if (!bKeepJacAcrossSteps) bUpdateJacobian = true;
     pNonlinearProblem = pNLP;
     VectorHandler* const pSol = pSolutionManager->pSolHdl();
     Size = pSol->iGetSize();
     pSol->Reset();
     DeltaX.ResizeReset(Size);
     XPrev.ResizeReset(Size);
     TmpRes.ResizeReset(Size);
}

void NoxNonlinearSolver::Residual(const VectorHandler* const pSol,
                                   VectorHandler* const pRes)
{
     CPUTimeGuard oCPUTimeRes(*this, CPU_RESIDUAL);
     // std::cerr << "DENTRO RESIDUAL" << std::endl;
     DeltaX.ScalarAddMul(*pSol, XPrev, -1.);
     XPrev = *pSol;
     pNonlinearProblem->Update(&DeltaX);
     pRes->Reset();
     VectorHandler* const pAbsRes = pGetResTest()->GetAbsRes();
     if (pAbsRes) pAbsRes->Reset();
     try {
          pNonlinearProblem->Residual(pRes, pAbsRes);
     } catch (const SolutionDataManager::ChangedEquationStructure&) {
          DEBUGCERR("Caught exception ChangedEquationStructure ...\n");
          if (bHonorJacRequest) {
               DEBUGCERR("Force update of preconditioner ...\n");
               ForcePrecondRebuild();
          }
     }
}

void NoxNonlinearSolver::Jacobian()
{
     CPUTimeGuard oCPUTimeJac(*this, CPU_JACOBIAN);
     bool bDone = false;
     MatrixHandler* pJac = nullptr;
     do {
          try {
               pJac = pSolutionManager->pMatHdl();
               pNonlinearProblem->Jacobian(pJac);
               pJac->PacMat();
               bDone = true;
          } catch (const MatrixHandler::ErrRebuildMatrix&) {
               silent_cout("NoxNonlinearSolver: rebuilding matrix...\n");
               pSolutionManager->MatrInitialize();
               // CrsMatrix pointer may have changed; re-wrap on next create_W_op()
               if (pModelEval) pModelEval->pJacobianOp = Teuchos::null;
          }
     } while (!bDone);
     ASSERT(pJac != nullptr);
#ifdef USE_MPI
     if (!bParallel || MBDynComm.Get_rank() == 0)
#endif
     {
          if (outputJac()) {
               silent_cout("Jacobian:\n");
               if (silent_out) pJac->Print(std::cout, MatrixHandler::MAT_PRINT_TRIPLET);
          }
     }
     ++TotJac;
}

void NoxNonlinearSolver::runPreLineSearch(const NOX::Solver::Generic&)
{
     DEBUGCERR("runPreLineSearch()\n");
     bInLineSearch = true;
     SetNonlinearSolverHint(LINESEARCH_LAMBDA_CURR, 1.);
     SetNonlinearSolverHint(LINESEARCH_ITERATION_CURR, 0);
}

void NoxNonlinearSolver::runPostLineSearch(const NOX::Solver::Generic&)
{
     DEBUGCERR("runPostLineSearch()\n");
     bInLineSearch = false;
     SetNonlinearSolverHint(LINESEARCH_LAMBDA_CURR, 1.);
     SetNonlinearSolverHint(LINESEARCH_ITERATION_CURR, 0);
}

} // anonymous namespace

/* =========================================================================
 * Public factory
 * ========================================================================= */
NonlinearSolver*
pAllocateNoxNonlinearSolver(const NonlinearSolverTestOptions& oSolverOpt,
                            const NoxSolverParameters& oParam)
{
     NoxNonlinearSolver* pNLS = nullptr;
     SAFENEWWITHCONSTRUCTOR(pNLS, NoxNonlinearSolver,
                            NoxNonlinearSolver(oSolverOpt, oParam));
     return pNLS;
}

#endif /* USE_TRILINOS */
