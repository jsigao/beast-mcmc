/*
 * CladeOperatedLoggerParser.java
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

package dr.evomodelxml.operators;

import java.io.PrintWriter;

import dr.evomodel.operators.AbstractTreeOperator;
import dr.evomodel.operators.CladeOperated;
import dr.evomodel.operators.CladeOperatedLogger;
import dr.inference.loggers.LogFormatter;
import dr.inference.loggers.TabDelimitedFormatter;
import dr.inferencexml.loggers.LoggerParser;
import dr.xml.*;

/**
 * @author Jiansi Gao
 */
public class CladeOperatedLoggerParser extends LoggerParser {
    public static final String LOG_CLADES_OPERATED = "logCladeOperated";
    public static final String APPEND = "append";

    public String getParserName() {
        return LOG_CLADES_OPERATED;
    }

    /**
     * @return an object based on the XML element it was passed.
     */
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        
        int logEvery = 0;
        if (xo.hasAttribute(LOG_EVERY)) {
            logEvery = xo.getIntegerAttribute(LOG_EVERY);
        }

        boolean append = false;
        if (xo.hasAttribute(APPEND)) {
            append = xo.getBooleanAttribute(APPEND);
        }

        String fileName = xo.getStringAttribute(FILE_NAME);

        final PrintWriter pw = XMLParser.getFilePrintWriter(xo, getParserName(), append);
        // final PrintWriter pw = getLogFile(xo, getParserName());
        final LogFormatter formatter = new TabDelimitedFormatter(pw);

        CladeOperated cladeOperated = new CladeOperated();

        for (int i = 0; i < xo.getChildCount(); i++) {
            Object child = xo.getChild(i);
            if (child instanceof AbstractTreeOperator) {
                ((AbstractTreeOperator) child).setCladeOperated(cladeOperated);
            }
        }

        return new CladeOperatedLogger(cladeOperated, formatter, logEvery, append, fileName);
    }

    public String getParserDescription() {
        return "Logs the clades that have been operated to a file.";
    }

    public Class getReturnType() {
        return CladeOperatedLogger.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }
    
    private final XMLSyntaxRule[] rules = {
            AttributeRule.newIntegerRule(LOG_EVERY, true),
            new StringAttributeRule(FILE_NAME,
                    "The name of the file to send log output to. " +
                            "If no file name is specified then log is sent to standard output", true),
            AttributeRule.newBooleanRule(APPEND, true),
            new ElementRule(AbstractTreeOperator.class, 1, Integer.MAX_VALUE)
    };
}
