/*
 * CladeOperatedLogger.java
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

package dr.evomodel.operators;

import dr.evomodel.operators.CladeOperated;

import dr.inference.loggers.LogFormatter;
import dr.inference.loggers.MCLogger;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.BitSet;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map;

/**
 * A logger for a random local tree variable.
 * It logs only the selected parameters in pairs of node number, variable value
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: TreeLogger.java,v 1.25 2006/09/05 13:29:34 rambaut Exp $
 */
public class CladeOperatedLogger extends MCLogger {

    private CladeOperated cladeOperated;
    private boolean append;
    private String fileName;

    public CladeOperatedLogger(CladeOperated cladeOperated, LogFormatter formatter, int logEvery, boolean append, String fileName) {
        super(formatter, logEvery, false);
        this.cladeOperated = cladeOperated;
        this.append = append;
        this.fileName = fileName;
    }

    public void startLogging() {
        boolean noHeader = false;
        File f = new File(fileName);
        if (f.exists() && append) {
            // System.out.println("clade appending");
            noHeader = readCladeFile(f);
        }
        if (noHeader == false) {
            String[] labels = {"state", "id", "clade"};
            logLabels(labels);
        }
    }

    public void log(long state) {

        if (logEvery <= 0 || ((state % logEvery) == 0)) {

            Map<BitSet, Integer> newCladeSet = cladeOperated.getNewCladeSet();
            String[] values = new String[3];
            values[0] = Long.toString(state);
            for (BitSet key : newCladeSet.keySet()) {
                values[1] = newCladeSet.get(key).toString();
                values[2] = key.toString();
                logValues(values);
            }

            cladeOperated.clearNewCladeSet();
        }
    }

    public BitSet parseBitSetString(String bitString) {
        
        int startIndex = bitString.indexOf('{');
        int endIndex = bitString.indexOf('}');
        if (endIndex == -1) {
            bitString = bitString.substring(startIndex + 1);
        } else {
            bitString = bitString.substring(startIndex + 1, endIndex);
        }

        String[] bitStringStr = bitString.split(", ");
        if (bitStringStr.length == 1) {
            bitStringStr = bitString.split(",");
        }

        BitSet bits = new BitSet();
        for (String str : bitStringStr) {
            int index = Integer.parseInt(str);
            bits.set(index);
        }

        return bits;
    }

    public boolean readCladeFile(File file) {
        
        try {
            BufferedReader reader = new BufferedReader(new FileReader(file));

            String line = reader.readLine();
            if (line == null) {
                reader.close();
                return false;
            }
            while ((line != null) && (line.startsWith("[") || line.startsWith("#"))) {
                line = reader.readLine();
            }
            
            if (line == null) {
                reader.close();
                return false;
            }
            String[] fields = line.split("\t");
            if (fields.length != 3 || fields[1].equals("id") == false || fields[2].equals("clade") == false) {
                reader.close();
                return false;
            }

            Map<BitSet, Integer> cladeSet = new HashMap<BitSet, Integer>();
            int cladeId = 0;
            line = reader.readLine();
            while (line != null) {
                fields = line.split("\t");
                if (fields.length != 3) {
                    if (cladeId == 0) {
                        reader.close();
                        return false;
                    } else {
                        throw new RuntimeException("Unable to read clade from the clade file");
                    }
                }
                // System.out.println("clade checking");
                cladeId = Integer.parseInt(fields[1]);
                BitSet bits = parseBitSetString(fields[2]);
                if (cladeSet.containsKey(bits)) {
                    // System.out.println("bits=" + bits.toString);
                    throw new RuntimeException("Duplicate clade found; Id=" + cladeId);
                } else {
                    cladeSet.put(bits, cladeId);
                }

                line = reader.readLine();
            }

            if (cladeSet.size() > 0) {
                cladeOperated.setCladeSet(cladeSet, cladeId);
            }

            reader.close();
            return true;

        } catch (IOException ioe) {
            throw new RuntimeException("Unable to read file: " + ioe.getMessage());
        }
    }

}
