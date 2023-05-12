/*
 * SiteLogLikelihoodLoggerParser.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.app.beagle.tools.parsers;

import dr.evomodel.treelikelihood.BeagleTreeLikelihood;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.inference.model.CompoundLikelihood;
import dr.app.beagle.tools.SiteLogLikelihoodLogger;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;
import dr.xml.XORRule;

public class SiteLogLikelihoodLoggerParser extends AbstractXMLObjectParser {

	public static final String SITE_LOGLIKELIHOOD_LOGGER = "siteLogLikelihood";

	@Override
	public Object parseXMLObject(XMLObject xo) throws XMLParseException {

		SiteLogLikelihoodLogger siteLogLikelihoodLogger;

		if (xo.getChild(BeagleTreeLikelihood.class) != null) {
			BeagleTreeLikelihood treeLikelihood = (BeagleTreeLikelihood) xo.getChild(BeagleTreeLikelihood.class);
			siteLogLikelihoodLogger = new SiteLogLikelihoodLogger(treeLikelihood);
		} else if (xo.getChild(TreeDataLikelihood.class) != null) {
			TreeDataLikelihood treeLikelihood = (TreeDataLikelihood) xo.getChild(TreeDataLikelihood.class);
			siteLogLikelihoodLogger = new SiteLogLikelihoodLogger(treeLikelihood);
		} else if (xo.getChild(CompoundLikelihood.class) != null) {
			CompoundLikelihood compoundLikelihood = (CompoundLikelihood) xo.getChild(CompoundLikelihood.class);
			siteLogLikelihoodLogger = new SiteLogLikelihoodLogger(compoundLikelihood);
		} else {
			throw new XMLParseException("No treeLikelihood or treeDatalikelihood available in siteLogLikelihood " + xo.getId());
		}

		return siteLogLikelihoodLogger;
	}// END: parseXMLObject

	@Override
	public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }// END: getSyntaxRules

	private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
		new XORRule(
			new XMLSyntaxRule[]{
				new ElementRule(BeagleTreeLikelihood.class),
				new ElementRule(TreeDataLikelihood.class),
				new ElementRule(CompoundLikelihood.class)
			}
		)
    };

	@Override
	public String getParserName() {
		return SITE_LOGLIKELIHOOD_LOGGER;
	}// END: getParserName

	@Override
	public String getParserDescription() {
		return "Beagle site logLikelihood";
	}// END: getParserDescription

	@Override
	public Class<SiteLogLikelihoodLogger> getReturnType() {
		return SiteLogLikelihoodLogger.class;
	}// getReturnType

}// END: class
