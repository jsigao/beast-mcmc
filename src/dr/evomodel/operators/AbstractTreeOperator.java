/*
 * AbstractTreeOperator.java
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

import dr.evomodel.tree.TreeModel;
import dr.evolution.tree.*;
import dr.inference.operators.SimpleMCMCOperator;

import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.loggers.NumberColumn;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public abstract class AbstractTreeOperator extends SimpleMCMCOperator {

	private long transitions = 0;

    protected CladeOperated cladeOperated;
    protected boolean logCladeOperated = false;

    protected long calculationCount = 0;

	/**
     * @return the number of transitions since last call to reset().
     */
    public long getTransitions() {
    	return transitions;
    }

    /**
     * Set the number of transitions since last call to reset(). This is used
     * to restore the state of the operator
     *
     * @param transitions number of transition
     */
    public void setTransitions(long transitions) {
    	this.transitions = transitions;
    }

    public void setCladeOperated(CladeOperated cladeOperated) {
    	this.cladeOperated = cladeOperated;
        logCladeOperated = true;
    }

    protected void setLogCladeOperated(boolean logCladeOperated) {
        this.logCladeOperated = logCladeOperated;
    }

    protected int getCladeIdx(Tree tree, NodeRef node) {
        return cladeOperated.getCladeIdx(tree, node);
    }

    public double getTransistionProbability() {
        final long accepted = getAcceptCount();
        final long rejected = getRejectCount();
        final long transition = getTransitions();
        return (double) transition / (double) (accepted + rejected);
    }

	/* exchange sub-trees whose root are i and j */
	protected void exchangeNodes(TreeModel tree, NodeRef i, NodeRef j,
	                             NodeRef iP, NodeRef jP) {

	    tree.beginTreeEdit();
	    tree.removeChild(iP, i);
	    tree.removeChild(jP, j);
	    tree.addChild(jP, i);
	    tree.addChild(iP, j);

        tree.endTreeEdit();
	}

	public void reset() {
        super.reset();
        transitions = 0;
    }

	/**
	 * @param tree   the tree
	 * @param parent the parent
	 * @param child  the child that you want the sister of
	 * @return the other child of the given parent.
	 */
    protected NodeRef getOtherChild(Tree tree, NodeRef parent, NodeRef child) {
        if( tree.getChild(parent, 0) == child ) {
            return tree.getChild(parent, 1);
        } else {
            return tree.getChild(parent, 0);
        }
    }

    protected int getNodeDistance(Tree tree, NodeRef i, NodeRef j) {
        int count = 0;

        while( i != j ) {
            count++;
            if( tree.getNodeHeight(i) < tree.getNodeHeight(j) ) {
                i = tree.getParent(i);
            } else {
                j = tree.getParent(j);
            }
        }
        return count;
    }

    public void addCalculationCount(long count) {
        super.addCalculationCount(count);
        calculationCount = count;
    }

    public String toStringStat(Integer[] v) {
		StringBuilder sb = new StringBuilder("{");
		for (int i = 0; i < v.length; i++) {
            sb.append(v[i]);
            if (i < v.length - 1) {
                sb.append(",");
            }
        }
        sb.append("}");
        return sb.toString();
	}

    public String toStringStat(Long[] v) {
		StringBuilder sb = new StringBuilder("{");
		for (int i = 0; i < v.length; i++) {
            sb.append(v[i]);
            if (i < v.length - 1) {
                sb.append(",");
            }
        }
        sb.append("}");
        return sb.toString();
	}

    public String toStringStat(Double[] v) {
		StringBuilder sb = new StringBuilder("{");
		for (int i = 0; i < v.length; i++) {
            NumberColumn column = new OperatorStatColumn(v[i].doubleValue());
            sb.append(column.getFormatted());
            if (i < v.length - 1) {
                sb.append(",");
            }
        }
        sb.append("}");
        return sb.toString();
	}

    protected class OperatorStatColumn extends NumberColumn {
        private final double value;

		public OperatorStatColumn(double value) {
			super("OperatorStatColumn", 4);
			this.value = value;
		}

		@Override
		public double getDoubleValue() {
			return value;
		}
	}

    public LogColumn getOperatorColumnInt(String label, List<Integer> statlist) {

        LogColumn column = new LogColumn.Abstract(getOperatorName() + "_" + label) {
            @Override
            protected String getFormattedValue() {
                Integer[] stats = statlist.toArray(new Integer[statlist.size()]);
                statlist.clear();
                return toStringStat(stats);
            }
        };

        return column;
    }

    public LogColumn getOperatorColumnLong(String label, List<Long> statlist) {

        LogColumn column = new LogColumn.Abstract(getOperatorName() + "_" + label) {
            @Override
            protected String getFormattedValue() {
                Long[] stats = statlist.toArray(new Long[statlist.size()]);
                statlist.clear();
                return toStringStat(stats);
            }
        };

        return column;
    }

    public LogColumn getOperatorColumnDouble(String label, List<Double> statlist) {

        LogColumn column = new LogColumn.Abstract(getOperatorName() + "_" + label) {
            @Override
            protected String getFormattedValue() {
                Double[] stats = statlist.toArray(new Double[statlist.size()]);
                statlist.clear();
                return toStringStat(stats);
            }
        };

        return column;
    }
}
