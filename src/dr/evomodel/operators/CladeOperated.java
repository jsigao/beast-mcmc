/*
 * CladeOperated.java
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

/**
 * 
 */
package dr.evomodel.operators;

import dr.evolution.tree.*;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * @author Jiansi Gao
 *
 */
public class CladeOperated {

    // private final List<AbstractTreeOperator> operators = new ArrayList<AbstractTreeOperator>();

    private Map<BitSet, Integer> cladeSet;
    private Map<BitSet, Integer> newCladeSet;
    private int cladeIdLast = 0;

    public CladeOperated() {
        cladeSet = new HashMap<BitSet, Integer>();
        newCladeSet = new LinkedHashMap<BitSet, Integer>();
    }

    // public void addOperator(AbstractTreeOperator op) {
	// 	// operators.add(op);
    //     op.setCladeOperated(this);
	// }

    protected BitSet getClade(Tree tree, NodeRef node) {

        // create a new bit set for this clade
        BitSet bits = new BitSet();

        // check if the node is external
        if (tree.isExternal(node)) {
            // if so, the only taxon in the clade am I
            String taxonId = tree.getNodeTaxon(node).getId();
            bits.set(tree.getTaxonIndex(taxonId));
            // int index = node.getNumber();
            // bits.set(index);
        } else {
            // otherwise, call all children and add its taxon together to one
            // clade
            for (int i = 0; i < tree.getChildCount(node); i++) {
                NodeRef child = tree.getChild(node, i);
                bits.or(getClade(tree, child));
            }
        }

        return bits;
    }

    public void clearNewCladeSet() {
        newCladeSet.clear();
    }

    public Map<BitSet, Integer> getNewCladeSet() {
        return newCladeSet;
    }

    public void setCladeSet(Map<BitSet, Integer> cladeSet, int cladeIdLast) {
        this.cladeSet = cladeSet;
        this.cladeIdLast = cladeIdLast;
    }

    public void updateCladeSet(BitSet bits) {

        if (!cladeSet.containsKey(bits)) {
            cladeIdLast++;
            cladeSet.put(bits, cladeIdLast);
            newCladeSet.put(bits, cladeIdLast);
        }
    }

    public int getCladeIdx(Tree tree, NodeRef node) {

        if (cladeIdLast == 0) {
            cladeIdLast = tree.getExternalNodeCount();
        }

        int index;
        if (tree.isExternal(node)) {
            String taxonId = tree.getNodeTaxon(node).getId();
            index = tree.getTaxonIndex(taxonId);
            // index = node.getNumber();
        } else {
            BitSet bits = getClade(tree, node);
            updateCladeSet(bits);
            index = cladeSet.get(bits);
        }

        return index;
    }

}
