/*
 * SiteLogLikelihoodLogger.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.app.beagle.tools;

import dr.evomodel.treelikelihood.BeagleTreeLikelihood;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.inference.model.CompoundLikelihood;
import dr.inference.model.Likelihood;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.loggers.NumberColumn;
import dr.xml.Reportable;

public class SiteLogLikelihoodLogger implements Loggable, Reportable {

	BeagleTreeLikelihood beagleTreeLikelihood;
	TreeDataLikelihood treeDataLikelihood;
	CompoundLikelihood compoundLikelihood;
	
	int patternCount;
	private SiteLogLikelihoodColumn[] columns = null;
	private double[] siteLikelihoods;
	private boolean useTreeLikelihood;
	private boolean useCompoundLikelihood = false;
	private int partitionCount = 1;

	public SiteLogLikelihoodLogger(BeagleTreeLikelihood beagleTreeLikelihood) {
		this.beagleTreeLikelihood = beagleTreeLikelihood;
		useTreeLikelihood = true;
	}// END: Constructor

	public SiteLogLikelihoodLogger(TreeDataLikelihood treeDataLikelihood) {
		this.treeDataLikelihood = treeDataLikelihood;
		useTreeLikelihood = false;
	}

	public SiteLogLikelihoodLogger(CompoundLikelihood compoundLikelihood) {
		this.compoundLikelihood = compoundLikelihood;
		useCompoundLikelihood = true;
		partitionCount = compoundLikelihood.getLikelihoodCount();
		if (compoundLikelihood.getLikelihood(0) instanceof BeagleTreeLikelihood) {
			useTreeLikelihood = true;
		} else {
			useTreeLikelihood = false;
		}
	}

	private void updateSiteLogLikelihood() {
		
		if (useCompoundLikelihood == false) {
			if (useTreeLikelihood) {
				siteLikelihoods = beagleTreeLikelihood.getSiteLogLikelihoods();
			} else {
				siteLikelihoods = treeDataLikelihood.getSiteLogLikelihoods();
			}
		} else {
			double[][] siteLikelihoodsTmp = new double[partitionCount][];
			int i = 0;
			int size = 0;
			int[] sizes = new int[partitionCount];
			for (Likelihood likelihood : compoundLikelihood.getLikelihoods()) {
				if (useTreeLikelihood) {
					siteLikelihoodsTmp[i] = ((BeagleTreeLikelihood) likelihood).getSiteLogLikelihoods();
				} else {
					siteLikelihoodsTmp[i] = ((TreeDataLikelihood) likelihood).getSiteLogLikelihoods();
				}
				sizes[i] = siteLikelihoodsTmp[i].length;
				size += sizes[i];
				i++;
			}
			
			siteLikelihoods = new double[size];
			size = 0;
			for (i = 0; i < partitionCount; i++) {
				System.arraycopy(siteLikelihoodsTmp[i], 0, siteLikelihoods, size, sizes[i]);
				size += sizes[i];
			}
		}
		
		patternCount = siteLikelihoods.length;
	}

	@Override
	public LogColumn[] getColumns() {

		updateSiteLogLikelihood();

		if (columns == null) {
			columns = new SiteLogLikelihoodColumn[patternCount];
			for (int site = 0; site < patternCount; site++) {
				columns[site] = new SiteLogLikelihoodColumn(site);
			}
		}

		return columns;
	}// END: getColumns

	private double getSiteLogLikelihood(int site) {
		if (site == 0) {
			updateSiteLogLikelihood();
		}
		return siteLikelihoods[site];
	}// END: getSiteLogLikelihoods

	private class SiteLogLikelihoodColumn extends NumberColumn {

		final int site;

		public SiteLogLikelihoodColumn(int site) {
			super("siteLikelihood" + "_" + (site + 1));
			this.site = site;
		}

		@Override
		public double getDoubleValue() {
			return getSiteLogLikelihood(site);
		}

	}// END: SiteLogLikelihoodColumn class

	public String toString() {
		LogColumn[] columns = getColumns();
		
		StringBuilder sb = new StringBuilder();
		for (int site = 0; site < patternCount; ++site) {
			if (site > 0) {
				sb.append(", ");
			}

			sb.append(columns[site].getFormatted());
		}

		return sb.toString();
	}// END: toString

	@Override
	public String getReport() {
		return toString();
	}// END: getReport

}// END: class
