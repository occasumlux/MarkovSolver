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
    mmmdialog.cpp \
    tableview.cpp \
    barchartwindow.cpp \
    qcustomplot.cpp \
    mer1dialog.cpp

HEADERS  += mainwindow.hpp \
    mmmdialog.hpp \
    ../eigen3/Eigen/src/Cholesky/LDLT.h \
    ../eigen3/Eigen/src/Cholesky/LLT.h \
    ../eigen3/Eigen/src/Cholesky/LLT_LAPACKE.h \
    ../eigen3/Eigen/src/CholmodSupport/CholmodSupport.h \
    ../eigen3/Eigen/src/Core/arch/AltiVec/Complex.h \
    ../eigen3/Eigen/src/Core/arch/AltiVec/MathFunctions.h \
    ../eigen3/Eigen/src/Core/arch/AltiVec/PacketMath.h \
    ../eigen3/Eigen/src/Core/arch/AVX/Complex.h \
    ../eigen3/Eigen/src/Core/arch/AVX/MathFunctions.h \
    ../eigen3/Eigen/src/Core/arch/AVX/PacketMath.h \
    ../eigen3/Eigen/src/Core/arch/AVX/TypeCasting.h \
    ../eigen3/Eigen/src/Core/arch/AVX512/MathFunctions.h \
    ../eigen3/Eigen/src/Core/arch/AVX512/PacketMath.h \
    ../eigen3/Eigen/src/Core/arch/CUDA/Complex.h \
    ../eigen3/Eigen/src/Core/arch/CUDA/Half.h \
    ../eigen3/Eigen/src/Core/arch/CUDA/MathFunctions.h \
    ../eigen3/Eigen/src/Core/arch/CUDA/PacketMath.h \
    ../eigen3/Eigen/src/Core/arch/CUDA/PacketMathHalf.h \
    ../eigen3/Eigen/src/Core/arch/CUDA/TypeCasting.h \
    ../eigen3/Eigen/src/Core/arch/Default/ConjHelper.h \
    ../eigen3/Eigen/src/Core/arch/Default/Settings.h \
    ../eigen3/Eigen/src/Core/arch/NEON/Complex.h \
    ../eigen3/Eigen/src/Core/arch/NEON/MathFunctions.h \
    ../eigen3/Eigen/src/Core/arch/NEON/PacketMath.h \
    ../eigen3/Eigen/src/Core/arch/SSE/Complex.h \
    ../eigen3/Eigen/src/Core/arch/SSE/MathFunctions.h \
    ../eigen3/Eigen/src/Core/arch/SSE/PacketMath.h \
    ../eigen3/Eigen/src/Core/arch/SSE/TypeCasting.h \
    ../eigen3/Eigen/src/Core/arch/ZVector/Complex.h \
    ../eigen3/Eigen/src/Core/arch/ZVector/MathFunctions.h \
    ../eigen3/Eigen/src/Core/arch/ZVector/PacketMath.h \
    ../eigen3/Eigen/src/Core/functors/AssignmentFunctors.h \
    ../eigen3/Eigen/src/Core/functors/BinaryFunctors.h \
    ../eigen3/Eigen/src/Core/functors/NullaryFunctors.h \
    ../eigen3/Eigen/src/Core/functors/StlFunctors.h \
    ../eigen3/Eigen/src/Core/functors/TernaryFunctors.h \
    ../eigen3/Eigen/src/Core/functors/UnaryFunctors.h \
    ../eigen3/Eigen/src/Core/products/GeneralBlockPanelKernel.h \
    ../eigen3/Eigen/src/Core/products/GeneralMatrixMatrix.h \
    ../eigen3/Eigen/src/Core/products/GeneralMatrixMatrix_BLAS.h \
    ../eigen3/Eigen/src/Core/products/GeneralMatrixMatrixTriangular.h \
    ../eigen3/Eigen/src/Core/products/GeneralMatrixMatrixTriangular_BLAS.h \
    ../eigen3/Eigen/src/Core/products/GeneralMatrixVector.h \
    ../eigen3/Eigen/src/Core/products/GeneralMatrixVector_BLAS.h \
    ../eigen3/Eigen/src/Core/products/Parallelizer.h \
    ../eigen3/Eigen/src/Core/products/SelfadjointMatrixMatrix.h \
    ../eigen3/Eigen/src/Core/products/SelfadjointMatrixMatrix_BLAS.h \
    ../eigen3/Eigen/src/Core/products/SelfadjointMatrixVector.h \
    ../eigen3/Eigen/src/Core/products/SelfadjointMatrixVector_BLAS.h \
    ../eigen3/Eigen/src/Core/products/SelfadjointProduct.h \
    ../eigen3/Eigen/src/Core/products/SelfadjointRank2Update.h \
    ../eigen3/Eigen/src/Core/products/TriangularMatrixMatrix.h \
    ../eigen3/Eigen/src/Core/products/TriangularMatrixMatrix_BLAS.h \
    ../eigen3/Eigen/src/Core/products/TriangularMatrixVector.h \
    ../eigen3/Eigen/src/Core/products/TriangularMatrixVector_BLAS.h \
    ../eigen3/Eigen/src/Core/products/TriangularSolverMatrix.h \
    ../eigen3/Eigen/src/Core/products/TriangularSolverMatrix_BLAS.h \
    ../eigen3/Eigen/src/Core/products/TriangularSolverVector.h \
    ../eigen3/Eigen/src/Core/util/BlasUtil.h \
    ../eigen3/Eigen/src/Core/util/Constants.h \
    ../eigen3/Eigen/src/Core/util/DisableStupidWarnings.h \
    ../eigen3/Eigen/src/Core/util/ForwardDeclarations.h \
    ../eigen3/Eigen/src/Core/util/Macros.h \
    ../eigen3/Eigen/src/Core/util/Memory.h \
    ../eigen3/Eigen/src/Core/util/Meta.h \
    ../eigen3/Eigen/src/Core/util/MKL_support.h \
    ../eigen3/Eigen/src/Core/util/NonMPL2.h \
    ../eigen3/Eigen/src/Core/util/ReenableStupidWarnings.h \
    ../eigen3/Eigen/src/Core/util/StaticAssert.h \
    ../eigen3/Eigen/src/Core/util/XprHelper.h \
    ../eigen3/Eigen/src/Core/Array.h \
    ../eigen3/Eigen/src/Core/ArrayBase.h \
    ../eigen3/Eigen/src/Core/ArrayWrapper.h \
    ../eigen3/Eigen/src/Core/Assign.h \
    ../eigen3/Eigen/src/Core/Assign_MKL.h \
    ../eigen3/Eigen/src/Core/AssignEvaluator.h \
    ../eigen3/Eigen/src/Core/BandMatrix.h \
    ../eigen3/Eigen/src/Core/Block.h \
    ../eigen3/Eigen/src/Core/BooleanRedux.h \
    ../eigen3/Eigen/src/Core/CommaInitializer.h \
    ../eigen3/Eigen/src/Core/ConditionEstimator.h \
    ../eigen3/Eigen/src/Core/CoreEvaluators.h \
    ../eigen3/Eigen/src/Core/CoreIterators.h \
    ../eigen3/Eigen/src/Core/CwiseBinaryOp.h \
    ../eigen3/Eigen/src/Core/CwiseNullaryOp.h \
    ../eigen3/Eigen/src/Core/CwiseTernaryOp.h \
    ../eigen3/Eigen/src/Core/CwiseUnaryOp.h \
    ../eigen3/Eigen/src/Core/CwiseUnaryView.h \
    ../eigen3/Eigen/src/Core/DenseBase.h \
    ../eigen3/Eigen/src/Core/DenseCoeffsBase.h \
    ../eigen3/Eigen/src/Core/DenseStorage.h \
    ../eigen3/Eigen/src/Core/Diagonal.h \
    ../eigen3/Eigen/src/Core/DiagonalMatrix.h \
    ../eigen3/Eigen/src/Core/DiagonalProduct.h \
    ../eigen3/Eigen/src/Core/Dot.h \
    ../eigen3/Eigen/src/Core/EigenBase.h \
    ../eigen3/Eigen/src/Core/ForceAlignedAccess.h \
    ../eigen3/Eigen/src/Core/Fuzzy.h \
    ../eigen3/Eigen/src/Core/GeneralProduct.h \
    ../eigen3/Eigen/src/Core/GenericPacketMath.h \
    ../eigen3/Eigen/src/Core/GlobalFunctions.h \
    ../eigen3/Eigen/src/Core/Inverse.h \
    ../eigen3/Eigen/src/Core/IO.h \
    ../eigen3/Eigen/src/Core/Map.h \
    ../eigen3/Eigen/src/Core/MapBase.h \
    ../eigen3/Eigen/src/Core/MathFunctions.h \
    ../eigen3/Eigen/src/Core/MathFunctionsImpl.h \
    ../eigen3/Eigen/src/Core/Matrix.h \
    ../eigen3/Eigen/src/Core/MatrixBase.h \
    ../eigen3/Eigen/src/Core/NestByValue.h \
    ../eigen3/Eigen/src/Core/NoAlias.h \
    ../eigen3/Eigen/src/Core/NumTraits.h \
    ../eigen3/Eigen/src/Core/PermutationMatrix.h \
    ../eigen3/Eigen/src/Core/PlainObjectBase.h \
    ../eigen3/Eigen/src/Core/Product.h \
    ../eigen3/Eigen/src/Core/ProductEvaluators.h \
    ../eigen3/Eigen/src/Core/Random.h \
    ../eigen3/Eigen/src/Core/Redux.h \
    ../eigen3/Eigen/src/Core/Ref.h \
    ../eigen3/Eigen/src/Core/Replicate.h \
    ../eigen3/Eigen/src/Core/ReturnByValue.h \
    ../eigen3/Eigen/src/Core/Reverse.h \
    ../eigen3/Eigen/src/Core/Select.h \
    ../eigen3/Eigen/src/Core/SelfAdjointView.h \
    ../eigen3/Eigen/src/Core/SelfCwiseBinaryOp.h \
    ../eigen3/Eigen/src/Core/Solve.h \
    ../eigen3/Eigen/src/Core/SolverBase.h \
    ../eigen3/Eigen/src/Core/SolveTriangular.h \
    ../eigen3/Eigen/src/Core/StableNorm.h \
    ../eigen3/Eigen/src/Core/Stride.h \
    ../eigen3/Eigen/src/Core/Swap.h \
    ../eigen3/Eigen/src/Core/Transpose.h \
    ../eigen3/Eigen/src/Core/Transpositions.h \
    ../eigen3/Eigen/src/Core/TriangularMatrix.h \
    ../eigen3/Eigen/src/Core/VectorBlock.h \
    ../eigen3/Eigen/src/Core/VectorwiseOp.h \
    ../eigen3/Eigen/src/Core/Visitor.h \
    ../eigen3/Eigen/src/Eigenvalues/ComplexEigenSolver.h \
    ../eigen3/Eigen/src/Eigenvalues/ComplexSchur.h \
    ../eigen3/Eigen/src/Eigenvalues/ComplexSchur_LAPACKE.h \
    ../eigen3/Eigen/src/Eigenvalues/EigenSolver.h \
    ../eigen3/Eigen/src/Eigenvalues/GeneralizedEigenSolver.h \
    ../eigen3/Eigen/src/Eigenvalues/GeneralizedSelfAdjointEigenSolver.h \
    ../eigen3/Eigen/src/Eigenvalues/HessenbergDecomposition.h \
    ../eigen3/Eigen/src/Eigenvalues/MatrixBaseEigenvalues.h \
    ../eigen3/Eigen/src/Eigenvalues/RealQZ.h \
    ../eigen3/Eigen/src/Eigenvalues/RealSchur.h \
    ../eigen3/Eigen/src/Eigenvalues/RealSchur_LAPACKE.h \
    ../eigen3/Eigen/src/Eigenvalues/SelfAdjointEigenSolver.h \
    ../eigen3/Eigen/src/Eigenvalues/SelfAdjointEigenSolver_LAPACKE.h \
    ../eigen3/Eigen/src/Eigenvalues/Tridiagonalization.h \
    ../eigen3/Eigen/src/Geometry/arch/Geometry_SSE.h \
    ../eigen3/Eigen/src/Geometry/AlignedBox.h \
    ../eigen3/Eigen/src/Geometry/AngleAxis.h \
    ../eigen3/Eigen/src/Geometry/EulerAngles.h \
    ../eigen3/Eigen/src/Geometry/Homogeneous.h \
    ../eigen3/Eigen/src/Geometry/Hyperplane.h \
    ../eigen3/Eigen/src/Geometry/OrthoMethods.h \
    ../eigen3/Eigen/src/Geometry/ParametrizedLine.h \
    ../eigen3/Eigen/src/Geometry/Quaternion.h \
    ../eigen3/Eigen/src/Geometry/Rotation2D.h \
    ../eigen3/Eigen/src/Geometry/RotationBase.h \
    ../eigen3/Eigen/src/Geometry/Scaling.h \
    ../eigen3/Eigen/src/Geometry/Transform.h \
    ../eigen3/Eigen/src/Geometry/Translation.h \
    ../eigen3/Eigen/src/Geometry/Umeyama.h \
    ../eigen3/Eigen/src/Householder/BlockHouseholder.h \
    ../eigen3/Eigen/src/Householder/Householder.h \
    ../eigen3/Eigen/src/Householder/HouseholderSequence.h \
    ../eigen3/Eigen/src/IterativeLinearSolvers/BasicPreconditioners.h \
    ../eigen3/Eigen/src/IterativeLinearSolvers/BiCGSTAB.h \
    ../eigen3/Eigen/src/IterativeLinearSolvers/ConjugateGradient.h \
    ../eigen3/Eigen/src/IterativeLinearSolvers/IncompleteCholesky.h \
    ../eigen3/Eigen/src/IterativeLinearSolvers/IncompleteLUT.h \
    ../eigen3/Eigen/src/IterativeLinearSolvers/IterativeSolverBase.h \
    ../eigen3/Eigen/src/IterativeLinearSolvers/LeastSquareConjugateGradient.h \
    ../eigen3/Eigen/src/IterativeLinearSolvers/SolveWithGuess.h \
    ../eigen3/Eigen/src/Jacobi/Jacobi.h \
    ../eigen3/Eigen/src/LU/arch/Inverse_SSE.h \
    ../eigen3/Eigen/src/LU/Determinant.h \
    ../eigen3/Eigen/src/LU/FullPivLU.h \
    ../eigen3/Eigen/src/LU/InverseImpl.h \
    ../eigen3/Eigen/src/LU/PartialPivLU.h \
    ../eigen3/Eigen/src/LU/PartialPivLU_LAPACKE.h \
    ../eigen3/Eigen/src/MetisSupport/MetisSupport.h \
    ../eigen3/Eigen/src/misc/blas.h \
    ../eigen3/Eigen/src/misc/Image.h \
    ../eigen3/Eigen/src/misc/Kernel.h \
    ../eigen3/Eigen/src/misc/lapack.h \
    ../eigen3/Eigen/src/misc/lapacke.h \
    ../eigen3/Eigen/src/misc/lapacke_mangling.h \
    ../eigen3/Eigen/src/misc/RealSvd2x2.h \
    ../eigen3/Eigen/src/OrderingMethods/Amd.h \
    ../eigen3/Eigen/src/OrderingMethods/Eigen_Colamd.h \
    ../eigen3/Eigen/src/OrderingMethods/Ordering.h \
    ../eigen3/Eigen/src/PardisoSupport/PardisoSupport.h \
    ../eigen3/Eigen/src/PaStiXSupport/PaStiXSupport.h \
    ../eigen3/Eigen/src/plugins/ArrayCwiseBinaryOps.h \
    ../eigen3/Eigen/src/plugins/ArrayCwiseUnaryOps.h \
    ../eigen3/Eigen/src/plugins/BlockMethods.h \
    ../eigen3/Eigen/src/plugins/CommonCwiseBinaryOps.h \
    ../eigen3/Eigen/src/plugins/CommonCwiseUnaryOps.h \
    ../eigen3/Eigen/src/plugins/MatrixCwiseBinaryOps.h \
    ../eigen3/Eigen/src/plugins/MatrixCwiseUnaryOps.h \
    ../eigen3/Eigen/src/QR/ColPivHouseholderQR.h \
    ../eigen3/Eigen/src/QR/ColPivHouseholderQR_LAPACKE.h \
    ../eigen3/Eigen/src/QR/CompleteOrthogonalDecomposition.h \
    ../eigen3/Eigen/src/QR/FullPivHouseholderQR.h \
    ../eigen3/Eigen/src/QR/HouseholderQR.h \
    ../eigen3/Eigen/src/QR/HouseholderQR_LAPACKE.h \
    ../eigen3/Eigen/src/SparseCholesky/SimplicialCholesky.h \
    ../eigen3/Eigen/src/SparseCholesky/SimplicialCholesky_impl.h \
    ../eigen3/Eigen/src/SparseCore/AmbiVector.h \
    ../eigen3/Eigen/src/SparseCore/CompressedStorage.h \
    ../eigen3/Eigen/src/SparseCore/ConservativeSparseSparseProduct.h \
    ../eigen3/Eigen/src/SparseCore/MappedSparseMatrix.h \
    ../eigen3/Eigen/src/SparseCore/SparseAssign.h \
    ../eigen3/Eigen/src/SparseCore/SparseBlock.h \
    ../eigen3/Eigen/src/SparseCore/SparseColEtree.h \
    ../eigen3/Eigen/src/SparseCore/SparseCompressedBase.h \
    ../eigen3/Eigen/src/SparseCore/SparseCwiseBinaryOp.h \
    ../eigen3/Eigen/src/SparseCore/SparseCwiseUnaryOp.h \
    ../eigen3/Eigen/src/SparseCore/SparseDenseProduct.h \
    ../eigen3/Eigen/src/SparseCore/SparseDiagonalProduct.h \
    ../eigen3/Eigen/src/SparseCore/SparseDot.h \
    ../eigen3/Eigen/src/SparseCore/SparseFuzzy.h \
    ../eigen3/Eigen/src/SparseCore/SparseMap.h \
    ../eigen3/Eigen/src/SparseCore/SparseMatrix.h \
    ../eigen3/Eigen/src/SparseCore/SparseMatrixBase.h \
    ../eigen3/Eigen/src/SparseCore/SparsePermutation.h \
    ../eigen3/Eigen/src/SparseCore/SparseProduct.h \
    ../eigen3/Eigen/src/SparseCore/SparseRedux.h \
    ../eigen3/Eigen/src/SparseCore/SparseRef.h \
    ../eigen3/Eigen/src/SparseCore/SparseSelfAdjointView.h \
    ../eigen3/Eigen/src/SparseCore/SparseSolverBase.h \
    ../eigen3/Eigen/src/SparseCore/SparseSparseProductWithPruning.h \
    ../eigen3/Eigen/src/SparseCore/SparseTranspose.h \
    ../eigen3/Eigen/src/SparseCore/SparseTriangularView.h \
    ../eigen3/Eigen/src/SparseCore/SparseUtil.h \
    ../eigen3/Eigen/src/SparseCore/SparseVector.h \
    ../eigen3/Eigen/src/SparseCore/SparseView.h \
    ../eigen3/Eigen/src/SparseCore/TriangularSolver.h \
    ../eigen3/Eigen/src/SparseLU/SparseLU.h \
    ../eigen3/Eigen/src/SparseLU/SparseLU_column_bmod.h \
    ../eigen3/Eigen/src/SparseLU/SparseLU_column_dfs.h \
    ../eigen3/Eigen/src/SparseLU/SparseLU_copy_to_ucol.h \
    ../eigen3/Eigen/src/SparseLU/SparseLU_gemm_kernel.h \
    ../eigen3/Eigen/src/SparseLU/SparseLU_heap_relax_snode.h \
    ../eigen3/Eigen/src/SparseLU/SparseLU_kernel_bmod.h \
    ../eigen3/Eigen/src/SparseLU/SparseLU_Memory.h \
    ../eigen3/Eigen/src/SparseLU/SparseLU_panel_bmod.h \
    ../eigen3/Eigen/src/SparseLU/SparseLU_panel_dfs.h \
    ../eigen3/Eigen/src/SparseLU/SparseLU_pivotL.h \
    ../eigen3/Eigen/src/SparseLU/SparseLU_pruneL.h \
    ../eigen3/Eigen/src/SparseLU/SparseLU_relax_snode.h \
    ../eigen3/Eigen/src/SparseLU/SparseLU_Structs.h \
    ../eigen3/Eigen/src/SparseLU/SparseLU_SupernodalMatrix.h \
    ../eigen3/Eigen/src/SparseLU/SparseLU_Utils.h \
    ../eigen3/Eigen/src/SparseLU/SparseLUImpl.h \
    ../eigen3/Eigen/src/SparseQR/SparseQR.h \
    ../eigen3/Eigen/src/SPQRSupport/SuiteSparseQRSupport.h \
    ../eigen3/Eigen/src/StlSupport/details.h \
    ../eigen3/Eigen/src/StlSupport/StdDeque.h \
    ../eigen3/Eigen/src/StlSupport/StdList.h \
    ../eigen3/Eigen/src/StlSupport/StdVector.h \
    ../eigen3/Eigen/src/SuperLUSupport/SuperLUSupport.h \
    ../eigen3/Eigen/src/SVD/BDCSVD.h \
    ../eigen3/Eigen/src/SVD/JacobiSVD.h \
    ../eigen3/Eigen/src/SVD/JacobiSVD_LAPACKE.h \
    ../eigen3/Eigen/src/SVD/SVDBase.h \
    ../eigen3/Eigen/src/SVD/UpperBidiagonalization.h \
    ../eigen3/Eigen/src/UmfPackSupport/UmfPackSupport.h \
    ../eigen3/Eigen/Cholesky \
    ../eigen3/Eigen/CholmodSupport \
    ../eigen3/Eigen/Core \
    ../eigen3/Eigen/Dense \
    ../eigen3/Eigen/Eigen \
    ../eigen3/Eigen/Eigenvalues \
    ../eigen3/Eigen/Geometry \
    ../eigen3/Eigen/Householder \
    ../eigen3/Eigen/IterativeLinearSolvers \
    ../eigen3/Eigen/Jacobi \
    ../eigen3/Eigen/LU \
    ../eigen3/Eigen/MetisSupport \
    ../eigen3/Eigen/OrderingMethods \
    ../eigen3/Eigen/PardisoSupport \
    ../eigen3/Eigen/PaStiXSupport \
    ../eigen3/Eigen/QR \
    ../eigen3/Eigen/QtAlignedMalloc \
    ../eigen3/Eigen/Sparse \
    ../eigen3/Eigen/SparseCholesky \
    ../eigen3/Eigen/SparseCore \
    ../eigen3/Eigen/SparseLU \
    ../eigen3/Eigen/SparseQR \
    ../eigen3/Eigen/SPQRSupport \
    ../eigen3/Eigen/StdDeque \
    ../eigen3/Eigen/StdList \
    ../eigen3/Eigen/StdVector \
    ../eigen3/Eigen/SuperLUSupport \
    ../eigen3/Eigen/SVD \
    ../eigen3/Eigen/UmfPackSupport \
    tableview.hpp \
    barchartwindow.hpp \
    qcustomplot.h \
    mer1dialog.hpp

FORMS    += mainwindow.ui \
    mmmdialog.ui \
    tableview.ui \
    barchartwindow.ui \
    mer1dialog.ui

DISTFILES += \
    ../eigen3/Eigen/CMakeLists.txt
