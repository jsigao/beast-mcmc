/*
 * GradientWrtParameterProvider.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.inference.hmc;

import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.preorder.NewTipFullConditionalDistributionDelegate;
import dr.evomodel.treedatalikelihood.preorder.NormalSufficientStatistics;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;

import java.util.List;

/**
 * @author Zhenyu Zhang
 * @author Marc A. Suchard
 */
public interface PrecisionVectorProductProvider {

    double[] getMultiplicationResultant(Parameter vector);

    class Generic implements PrecisionVectorProductProvider {

        private final MatrixParameterInterface matrix;

        public Generic(MatrixParameterInterface matrix) {
            this.matrix = matrix;
        }

        @Override
        public double[] getMultiplicationResultant(Parameter vector) {

            final int nRows = matrix.getRowDimension();
            final int nCols = matrix.getColumnDimension();

            assert (vector.getDimension() == nCols);

            double[] result = new double[nRows];

            for (int row = 0; row < nRows; ++row) {
                double sum = 0.0;
                for (int col = 0; col < nCols; ++col) {
                    sum += matrix.getParameterValue(row, col) * vector.getParameterValue(col);
                }
                result[row] = sum;
            }

            return result;
        }
    }

    // TODO Depends on evomodel objects, should move into evomodel package (limit cyclic dependencies)

}