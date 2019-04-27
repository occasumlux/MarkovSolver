#-------------------------------------------------
#
# Project created by QtCreator 2019-03-26T12:41:18
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets printsupport

TARGET = Markov
TEMPLATE = app

# The following define makes your compiler emit warnings if you use
# any feature of Qt which as been marked as deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0


SOURCES += main.cpp\
        mainwindow.cpp \
    tableview.cpp \
    barchartwindow.cpp \
    qcustomplot.cpp

HEADERS  += mainwindow.hpp \
    ../eigen-eigen-323c052e1731/Eigen/src/Cholesky/LDLT.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Cholesky/LLT.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Cholesky/LLT_LAPACKE.h \
    ../eigen-eigen-323c052e1731/Eigen/src/CholmodSupport/CholmodSupport.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/arch/AltiVec/Complex.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/arch/AltiVec/MathFunctions.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/arch/AltiVec/PacketMath.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/arch/AVX/Complex.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/arch/AVX/MathFunctions.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/arch/AVX/PacketMath.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/arch/AVX/TypeCasting.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/arch/AVX512/MathFunctions.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/arch/AVX512/PacketMath.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/arch/CUDA/Complex.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/arch/CUDA/Half.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/arch/CUDA/MathFunctions.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/arch/CUDA/PacketMath.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/arch/CUDA/PacketMathHalf.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/arch/CUDA/TypeCasting.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/arch/Default/ConjHelper.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/arch/Default/Settings.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/arch/NEON/Complex.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/arch/NEON/MathFunctions.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/arch/NEON/PacketMath.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/arch/SSE/Complex.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/arch/SSE/MathFunctions.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/arch/SSE/PacketMath.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/arch/SSE/TypeCasting.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/arch/ZVector/Complex.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/arch/ZVector/MathFunctions.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/arch/ZVector/PacketMath.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/functors/AssignmentFunctors.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/functors/BinaryFunctors.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/functors/NullaryFunctors.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/functors/StlFunctors.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/functors/TernaryFunctors.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/functors/UnaryFunctors.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/products/GeneralBlockPanelKernel.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/products/GeneralMatrixMatrix.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/products/GeneralMatrixMatrix_BLAS.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/products/GeneralMatrixMatrixTriangular.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/products/GeneralMatrixMatrixTriangular_BLAS.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/products/GeneralMatrixVector.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/products/GeneralMatrixVector_BLAS.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/products/Parallelizer.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/products/SelfadjointMatrixMatrix.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/products/SelfadjointMatrixMatrix_BLAS.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/products/SelfadjointMatrixVector.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/products/SelfadjointMatrixVector_BLAS.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/products/SelfadjointProduct.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/products/SelfadjointRank2Update.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/products/TriangularMatrixMatrix.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/products/TriangularMatrixMatrix_BLAS.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/products/TriangularMatrixVector.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/products/TriangularMatrixVector_BLAS.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/products/TriangularSolverMatrix.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/products/TriangularSolverMatrix_BLAS.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/products/TriangularSolverVector.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/util/BlasUtil.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/util/Constants.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/util/DisableStupidWarnings.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/util/ForwardDeclarations.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/util/Macros.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/util/Memory.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/util/Meta.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/util/MKL_support.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/util/NonMPL2.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/util/ReenableStupidWarnings.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/util/StaticAssert.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/util/XprHelper.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/Array.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/ArrayBase.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/ArrayWrapper.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/Assign.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/Assign_MKL.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/AssignEvaluator.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/BandMatrix.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/Block.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/BooleanRedux.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/CommaInitializer.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/ConditionEstimator.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/CoreEvaluators.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/CoreIterators.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/CwiseBinaryOp.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/CwiseNullaryOp.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/CwiseTernaryOp.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/CwiseUnaryOp.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/CwiseUnaryView.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/DenseBase.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/DenseCoeffsBase.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/DenseStorage.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/Diagonal.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/DiagonalMatrix.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/DiagonalProduct.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/Dot.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/EigenBase.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/ForceAlignedAccess.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/Fuzzy.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/GeneralProduct.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/GenericPacketMath.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/GlobalFunctions.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/Inverse.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/IO.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/Map.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/MapBase.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/MathFunctions.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/MathFunctionsImpl.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/Matrix.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/MatrixBase.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/NestByValue.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/NoAlias.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/NumTraits.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/PermutationMatrix.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/PlainObjectBase.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/Product.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/ProductEvaluators.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/Random.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/Redux.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/Ref.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/Replicate.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/ReturnByValue.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/Reverse.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/Select.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/SelfAdjointView.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/SelfCwiseBinaryOp.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/Solve.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/SolverBase.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/SolveTriangular.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/StableNorm.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/Stride.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/Swap.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/Transpose.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/Transpositions.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/TriangularMatrix.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/VectorBlock.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/VectorwiseOp.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Core/Visitor.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Eigenvalues/ComplexEigenSolver.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Eigenvalues/ComplexSchur.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Eigenvalues/ComplexSchur_LAPACKE.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Eigenvalues/EigenSolver.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Eigenvalues/GeneralizedEigenSolver.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Eigenvalues/GeneralizedSelfAdjointEigenSolver.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Eigenvalues/HessenbergDecomposition.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Eigenvalues/MatrixBaseEigenvalues.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Eigenvalues/RealQZ.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Eigenvalues/RealSchur.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Eigenvalues/RealSchur_LAPACKE.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Eigenvalues/SelfAdjointEigenSolver.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Eigenvalues/SelfAdjointEigenSolver_LAPACKE.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Eigenvalues/Tridiagonalization.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Geometry/arch/Geometry_SSE.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Geometry/AlignedBox.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Geometry/AngleAxis.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Geometry/EulerAngles.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Geometry/Homogeneous.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Geometry/Hyperplane.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Geometry/OrthoMethods.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Geometry/ParametrizedLine.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Geometry/Quaternion.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Geometry/Rotation2D.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Geometry/RotationBase.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Geometry/Scaling.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Geometry/Transform.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Geometry/Translation.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Geometry/Umeyama.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Householder/BlockHouseholder.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Householder/Householder.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Householder/HouseholderSequence.h \
    ../eigen-eigen-323c052e1731/Eigen/src/IterativeLinearSolvers/BasicPreconditioners.h \
    ../eigen-eigen-323c052e1731/Eigen/src/IterativeLinearSolvers/BiCGSTAB.h \
    ../eigen-eigen-323c052e1731/Eigen/src/IterativeLinearSolvers/ConjugateGradient.h \
    ../eigen-eigen-323c052e1731/Eigen/src/IterativeLinearSolvers/IncompleteCholesky.h \
    ../eigen-eigen-323c052e1731/Eigen/src/IterativeLinearSolvers/IncompleteLUT.h \
    ../eigen-eigen-323c052e1731/Eigen/src/IterativeLinearSolvers/IterativeSolverBase.h \
    ../eigen-eigen-323c052e1731/Eigen/src/IterativeLinearSolvers/LeastSquareConjugateGradient.h \
    ../eigen-eigen-323c052e1731/Eigen/src/IterativeLinearSolvers/SolveWithGuess.h \
    ../eigen-eigen-323c052e1731/Eigen/src/Jacobi/Jacobi.h \
    ../eigen-eigen-323c052e1731/Eigen/src/LU/arch/Inverse_SSE.h \
    ../eigen-eigen-323c052e1731/Eigen/src/LU/Determinant.h \
    ../eigen-eigen-323c052e1731/Eigen/src/LU/FullPivLU.h \
    ../eigen-eigen-323c052e1731/Eigen/src/LU/InverseImpl.h \
    ../eigen-eigen-323c052e1731/Eigen/src/LU/PartialPivLU.h \
    ../eigen-eigen-323c052e1731/Eigen/src/LU/PartialPivLU_LAPACKE.h \
    ../eigen-eigen-323c052e1731/Eigen/src/MetisSupport/MetisSupport.h \
    ../eigen-eigen-323c052e1731/Eigen/src/misc/blas.h \
    ../eigen-eigen-323c052e1731/Eigen/src/misc/Image.h \
    ../eigen-eigen-323c052e1731/Eigen/src/misc/Kernel.h \
    ../eigen-eigen-323c052e1731/Eigen/src/misc/lapack.h \
    ../eigen-eigen-323c052e1731/Eigen/src/misc/lapacke.h \
    ../eigen-eigen-323c052e1731/Eigen/src/misc/lapacke_mangling.h \
    ../eigen-eigen-323c052e1731/Eigen/src/misc/RealSvd2x2.h \
    ../eigen-eigen-323c052e1731/Eigen/src/OrderingMethods/Amd.h \
    ../eigen-eigen-323c052e1731/Eigen/src/OrderingMethods/Eigen_Colamd.h \
    ../eigen-eigen-323c052e1731/Eigen/src/OrderingMethods/Ordering.h \
    ../eigen-eigen-323c052e1731/Eigen/src/PardisoSupport/PardisoSupport.h \
    ../eigen-eigen-323c052e1731/Eigen/src/PaStiXSupport/PaStiXSupport.h \
    ../eigen-eigen-323c052e1731/Eigen/src/plugins/ArrayCwiseBinaryOps.h \
    ../eigen-eigen-323c052e1731/Eigen/src/plugins/ArrayCwiseUnaryOps.h \
    ../eigen-eigen-323c052e1731/Eigen/src/plugins/BlockMethods.h \
    ../eigen-eigen-323c052e1731/Eigen/src/plugins/CommonCwiseBinaryOps.h \
    ../eigen-eigen-323c052e1731/Eigen/src/plugins/CommonCwiseUnaryOps.h \
    ../eigen-eigen-323c052e1731/Eigen/src/plugins/MatrixCwiseBinaryOps.h \
    ../eigen-eigen-323c052e1731/Eigen/src/plugins/MatrixCwiseUnaryOps.h \
    ../eigen-eigen-323c052e1731/Eigen/src/QR/ColPivHouseholderQR.h \
    ../eigen-eigen-323c052e1731/Eigen/src/QR/ColPivHouseholderQR_LAPACKE.h \
    ../eigen-eigen-323c052e1731/Eigen/src/QR/CompleteOrthogonalDecomposition.h \
    ../eigen-eigen-323c052e1731/Eigen/src/QR/FullPivHouseholderQR.h \
    ../eigen-eigen-323c052e1731/Eigen/src/QR/HouseholderQR.h \
    ../eigen-eigen-323c052e1731/Eigen/src/QR/HouseholderQR_LAPACKE.h \
    ../eigen-eigen-323c052e1731/Eigen/src/SparseCholesky/SimplicialCholesky.h \
    ../eigen-eigen-323c052e1731/Eigen/src/SparseCholesky/SimplicialCholesky_impl.h \
    ../eigen-eigen-323c052e1731/Eigen/src/SparseCore/AmbiVector.h \
    ../eigen-eigen-323c052e1731/Eigen/src/SparseCore/CompressedStorage.h \
    ../eigen-eigen-323c052e1731/Eigen/src/SparseCore/ConservativeSparseSparseProduct.h \
    ../eigen-eigen-323c052e1731/Eigen/src/SparseCore/MappedSparseMatrix.h \
    ../eigen-eigen-323c052e1731/Eigen/src/SparseCore/SparseAssign.h \
    ../eigen-eigen-323c052e1731/Eigen/src/SparseCore/SparseBlock.h \
    ../eigen-eigen-323c052e1731/Eigen/src/SparseCore/SparseColEtree.h \
    ../eigen-eigen-323c052e1731/Eigen/src/SparseCore/SparseCompressedBase.h \
    ../eigen-eigen-323c052e1731/Eigen/src/SparseCore/SparseCwiseBinaryOp.h \
    ../eigen-eigen-323c052e1731/Eigen/src/SparseCore/SparseCwiseUnaryOp.h \
    ../eigen-eigen-323c052e1731/Eigen/src/SparseCore/SparseDenseProduct.h \
    ../eigen-eigen-323c052e1731/Eigen/src/SparseCore/SparseDiagonalProduct.h \
    ../eigen-eigen-323c052e1731/Eigen/src/SparseCore/SparseDot.h \
    ../eigen-eigen-323c052e1731/Eigen/src/SparseCore/SparseFuzzy.h \
    ../eigen-eigen-323c052e1731/Eigen/src/SparseCore/SparseMap.h \
    ../eigen-eigen-323c052e1731/Eigen/src/SparseCore/SparseMatrix.h \
    ../eigen-eigen-323c052e1731/Eigen/src/SparseCore/SparseMatrixBase.h \
    ../eigen-eigen-323c052e1731/Eigen/src/SparseCore/SparsePermutation.h \
    ../eigen-eigen-323c052e1731/Eigen/src/SparseCore/SparseProduct.h \
    ../eigen-eigen-323c052e1731/Eigen/src/SparseCore/SparseRedux.h \
    ../eigen-eigen-323c052e1731/Eigen/src/SparseCore/SparseRef.h \
    ../eigen-eigen-323c052e1731/Eigen/src/SparseCore/SparseSelfAdjointView.h \
    ../eigen-eigen-323c052e1731/Eigen/src/SparseCore/SparseSolverBase.h \
    ../eigen-eigen-323c052e1731/Eigen/src/SparseCore/SparseSparseProductWithPruning.h \
    ../eigen-eigen-323c052e1731/Eigen/src/SparseCore/SparseTranspose.h \
    ../eigen-eigen-323c052e1731/Eigen/src/SparseCore/SparseTriangularView.h \
    ../eigen-eigen-323c052e1731/Eigen/src/SparseCore/SparseUtil.h \
    ../eigen-eigen-323c052e1731/Eigen/src/SparseCore/SparseVector.h \
    ../eigen-eigen-323c052e1731/Eigen/src/SparseCore/SparseView.h \
    ../eigen-eigen-323c052e1731/Eigen/src/SparseCore/TriangularSolver.h \
    ../eigen-eigen-323c052e1731/Eigen/src/SparseLU/SparseLU.h \
    ../eigen-eigen-323c052e1731/Eigen/src/SparseLU/SparseLU_column_bmod.h \
    ../eigen-eigen-323c052e1731/Eigen/src/SparseLU/SparseLU_column_dfs.h \
    ../eigen-eigen-323c052e1731/Eigen/src/SparseLU/SparseLU_copy_to_ucol.h \
    ../eigen-eigen-323c052e1731/Eigen/src/SparseLU/SparseLU_gemm_kernel.h \
    ../eigen-eigen-323c052e1731/Eigen/src/SparseLU/SparseLU_heap_relax_snode.h \
    ../eigen-eigen-323c052e1731/Eigen/src/SparseLU/SparseLU_kernel_bmod.h \
    ../eigen-eigen-323c052e1731/Eigen/src/SparseLU/SparseLU_Memory.h \
    ../eigen-eigen-323c052e1731/Eigen/src/SparseLU/SparseLU_panel_bmod.h \
    ../eigen-eigen-323c052e1731/Eigen/src/SparseLU/SparseLU_panel_dfs.h \
    ../eigen-eigen-323c052e1731/Eigen/src/SparseLU/SparseLU_pivotL.h \
    ../eigen-eigen-323c052e1731/Eigen/src/SparseLU/SparseLU_pruneL.h \
    ../eigen-eigen-323c052e1731/Eigen/src/SparseLU/SparseLU_relax_snode.h \
    ../eigen-eigen-323c052e1731/Eigen/src/SparseLU/SparseLU_Structs.h \
    ../eigen-eigen-323c052e1731/Eigen/src/SparseLU/SparseLU_SupernodalMatrix.h \
    ../eigen-eigen-323c052e1731/Eigen/src/SparseLU/SparseLU_Utils.h \
    ../eigen-eigen-323c052e1731/Eigen/src/SparseLU/SparseLUImpl.h \
    ../eigen-eigen-323c052e1731/Eigen/src/SparseQR/SparseQR.h \
    ../eigen-eigen-323c052e1731/Eigen/src/SPQRSupport/SuiteSparseQRSupport.h \
    ../eigen-eigen-323c052e1731/Eigen/src/StlSupport/details.h \
    ../eigen-eigen-323c052e1731/Eigen/src/StlSupport/StdDeque.h \
    ../eigen-eigen-323c052e1731/Eigen/src/StlSupport/StdList.h \
    ../eigen-eigen-323c052e1731/Eigen/src/StlSupport/StdVector.h \
    ../eigen-eigen-323c052e1731/Eigen/src/SuperLUSupport/SuperLUSupport.h \
    ../eigen-eigen-323c052e1731/Eigen/src/SVD/BDCSVD.h \
    ../eigen-eigen-323c052e1731/Eigen/src/SVD/JacobiSVD.h \
    ../eigen-eigen-323c052e1731/Eigen/src/SVD/JacobiSVD_LAPACKE.h \
    ../eigen-eigen-323c052e1731/Eigen/src/SVD/SVDBase.h \
    ../eigen-eigen-323c052e1731/Eigen/src/SVD/UpperBidiagonalization.h \
    ../eigen-eigen-323c052e1731/Eigen/src/UmfPackSupport/UmfPackSupport.h \
    ../eigen-eigen-323c052e1731/Eigen/Cholesky \
    ../eigen-eigen-323c052e1731/Eigen/CholmodSupport \
    ../eigen-eigen-323c052e1731/Eigen/Core \
    ../eigen-eigen-323c052e1731/Eigen/Dense \
    ../eigen-eigen-323c052e1731/Eigen/Eigen \
    ../eigen-eigen-323c052e1731/Eigen/Eigenvalues \
    ../eigen-eigen-323c052e1731/Eigen/Geometry \
    ../eigen-eigen-323c052e1731/Eigen/Householder \
    ../eigen-eigen-323c052e1731/Eigen/IterativeLinearSolvers \
    ../eigen-eigen-323c052e1731/Eigen/Jacobi \
    ../eigen-eigen-323c052e1731/Eigen/LU \
    ../eigen-eigen-323c052e1731/Eigen/MetisSupport \
    ../eigen-eigen-323c052e1731/Eigen/OrderingMethods \
    ../eigen-eigen-323c052e1731/Eigen/PardisoSupport \
    ../eigen-eigen-323c052e1731/Eigen/PaStiXSupport \
    ../eigen-eigen-323c052e1731/Eigen/QR \
    ../eigen-eigen-323c052e1731/Eigen/QtAlignedMalloc \
    ../eigen-eigen-323c052e1731/Eigen/Sparse \
    ../eigen-eigen-323c052e1731/Eigen/SparseCholesky \
    ../eigen-eigen-323c052e1731/Eigen/SparseCore \
    ../eigen-eigen-323c052e1731/Eigen/SparseLU \
    ../eigen-eigen-323c052e1731/Eigen/SparseQR \
    ../eigen-eigen-323c052e1731/Eigen/SPQRSupport \
    ../eigen-eigen-323c052e1731/Eigen/StdDeque \
    ../eigen-eigen-323c052e1731/Eigen/StdList \
    ../eigen-eigen-323c052e1731/Eigen/StdVector \
    ../eigen-eigen-323c052e1731/Eigen/SuperLUSupport \
    ../eigen-eigen-323c052e1731/Eigen/SVD \
    ../eigen-eigen-323c052e1731/Eigen/UmfPackSupport \
    tableview.hpp \
    barchartwindow.hpp \
    qcustomplot.h

FORMS    += mainwindow.ui \
    tableview.ui \
    barchartwindow.ui

DISTFILES += \
    ../eigen-eigen-323c052e1731/Eigen/CMakeLists.txt
