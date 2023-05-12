/*
 * ExchangeOperator.java
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

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.TreeUtils;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treelikelihood.thorneytreelikelihood.ConstrainableTreeOperator;
import dr.math.MathUtils;

import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;

import java.util.ArrayDeque;
import java.util.Deque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.logging.Logger;

/**
 * Implements branch exchange operations. There is a NARROW and WIDE variety.
 * The narrow exchange is very similar to a rooted-tree nearest-neighbour
 * interchange but with the restriction that node height must remain consistent.
 * <p/>
 * KNOWN BUGS: WIDE operator cannot be used on trees with 4 or fewer tips!
 */
public class ExchangeOperator extends AbstractTreeOperator implements ConstrainableTreeOperator {

    public static final int NARROW = 0;
    public static final int WIDE = 1;

    private static final int MAX_TRIES = 100;

    private int mode = NARROW;
    private final TreeModel tree;

    private double[] distances;

    private int nodeDistance = 3;
    private double pathLength;
    private double[] nodeHeights = new double[2];
    private int[] leafCounts = new int[2];

    private List<Integer> nodeDistanceAccept;
    private List<Integer> nodeDistanceReject;
    private List<Double> pathLengthAccept;
    private List<Double> pathLengthReject;
    private List<Double> nodeHeight0Accept;
    private List<Double> nodeHeight0Reject;
    private List<Double> nodeHeight1Accept;
    private List<Double> nodeHeight1Reject;
    private List<Integer> leafCounts0Accept;
    private List<Integer> leafCounts0Reject;
    private List<Integer> leafCounts1Accept;
    private List<Integer> leafCounts1Reject;

    public ExchangeOperator(int mode, TreeModel tree, double weight) {
        this.mode = mode;
        this.tree = tree;
        setWeight(weight);

        if (mode != NARROW) {
            nodeDistanceAccept = new ArrayList<Integer>();
            nodeDistanceReject = new ArrayList<Integer>();
        }
        pathLengthAccept = new ArrayList<Double>();
        pathLengthReject = new ArrayList<Double>();
        nodeHeight0Accept = new ArrayList<Double>();
        nodeHeight0Reject = new ArrayList<Double>();
        nodeHeight1Accept = new ArrayList<Double>();
        nodeHeight1Reject = new ArrayList<Double>();
        leafCounts0Accept = new ArrayList<Integer>();
        leafCounts0Reject = new ArrayList<Integer>();
        leafCounts1Accept = new ArrayList<Integer>();
        leafCounts1Reject = new ArrayList<Integer>();
    }


    public double doOperation(TreeModel tree) {

        final int tipCount = tree.getExternalNodeCount();

        double hastingsRatio;

        switch( mode ) {
            case NARROW:
                hastingsRatio = (narrow(tree) ? 0.0 : Double.NEGATIVE_INFINITY);
                break;
            case WIDE:
                hastingsRatio = (wide(tree) ? 0.0 : Double.NEGATIVE_INFINITY);
                break;
            default:
                throw new IllegalArgumentException("Unknow Exchange Mode");
        }

        assert tree.getExternalNodeCount() == tipCount :
                "Lost some tips in " + ((mode == NARROW) ? "NARROW mode." : "WIDE mode.");

        return hastingsRatio;
    }
    public double doOperation(){
        return doOperation(tree);
    }

    /**
     * WARNING: Assumes strictly bifurcating tree.
     */
    public boolean narrow(TreeModel tree) {
        final int nNodes = tree.getNodeCount();
        final NodeRef root = tree.getRoot();

        NodeRef i = root;

        while( root == i || tree.getParent(i) == root ) {
            i = tree.getNode(MathUtils.nextInt(nNodes));
        }

        final NodeRef iParent = tree.getParent(i);
        final NodeRef iGrandParent = tree.getParent(iParent);
        NodeRef iUncle = tree.getChild(iGrandParent, 0);
        if( iUncle == iParent ) {
            iUncle = tree.getChild(iGrandParent, 1);
        }
        assert iUncle == getOtherChild(tree, iGrandParent, iParent);

        assert tree.getNodeHeight(i) <= tree.getNodeHeight(iGrandParent);

        pathLength = tree.getBranchLength(i) + tree.getBranchLength(iParent) + tree.getBranchLength(iUncle);
        nodeHeights[0] = tree.getNodeHeight(i);
        nodeHeights[1] = tree.getNodeHeight(iUncle);
        leafCounts[0] = TreeUtils.getLeafCount(tree, i);
        leafCounts[1] = TreeUtils.getLeafCount(tree, iUncle);

        if( tree.getNodeHeight(iUncle) < tree.getNodeHeight(iParent) ) {
            exchangeNodes(tree, i, iUncle, iParent, iGrandParent);

            // exchangeNodes generates the events
            //tree.pushTreeChangedEvent(iParent);
            //tree.pushTreeChangedEvent(iGrandParent);
            return true;
        }

        return false;
    }

    /**
     * WARNING: Assumes strictly bifurcating tree.
     */
    public boolean wide(TreeModel tree) {

        final int nodeCount = tree.getNodeCount();
        final NodeRef root = tree.getRoot();

        NodeRef i = root;

        while( root == i ) {
            i = tree.getNode(MathUtils.nextInt(nodeCount));
        }

        NodeRef j = i;
        while( j == i || j == root ) {
            j = tree.getNode(MathUtils.nextInt(nodeCount));
        }

        final NodeRef iP = tree.getParent(i);
        final NodeRef jP = tree.getParent(j);

        nodeDistance = getNodeDistance(tree, i, j);
        pathLength = TreeUtils.getPathLength(tree, i, j);
        nodeHeights[0] = tree.getNodeHeight(i);
        nodeHeights[1] = tree.getNodeHeight(j);
        leafCounts[0] = TreeUtils.getLeafCount(tree, i);
        leafCounts[1] = TreeUtils.getLeafCount(tree, j);

        if( (iP != jP) && (i != jP) && (j != iP)
                && (tree.getNodeHeight(j) < tree.getNodeHeight(iP))
                && (tree.getNodeHeight(i) < tree.getNodeHeight(jP)) ) {
            exchangeNodes(tree, i, j, iP, jP);
            // System.out.println("tries = " + tries+1);
            return true;
        }

        return false;
    }

    /* why not use Arrays.asList(a).indexOf(n) ? */
    private int indexOf(NodeRef[] a, NodeRef n) {

        for(int i = 0; i < a.length; i++) {
            if( a[i] == n ) {
                return i;
            }
        }
        return -1;
    }

    private double getWinningChance(int index) {

        double sum = 0;
        for( double distance : distances ) {
            sum += (1.0 / distance);
        }

        return (1.0 / distances[index]) / sum;

    }

    private void calcDistances(NodeRef[] nodes, NodeRef ref) {
        distances = new double[nodes.length];
        for(int i = 0; i < nodes.length; i++) {
            distances[i] = getNodeDistance(tree, ref, nodes[i]) + 1;
        }
    }

    private NodeRef getRandomNode(NodeRef[] nodes, NodeRef ref) {

        calcDistances(nodes, ref);
        double sum = 0;
        for( double distance : distances ) {
            sum += 1.0 / distance;
        }

        double randomValue = MathUtils.nextDouble() * sum;
        NodeRef n = null;
        for(int i = 0; i < distances.length; i++) {
            randomValue -= 1.0 / distances[i];

            if( randomValue <= 0 ) {
                n = nodes[i];
                break;
            }
        }
        return n;
    }

    public int getMode() {
        return mode;
    }

    public String getOperatorName() {
        return ((mode == NARROW) ? "Narrow" : "Wide") + " Exchange" + "(" + tree.getId() + ")";
    }

    public void accept(double deviation) {
        super.accept(deviation);

        if (mode != NARROW) {
            nodeDistanceAccept.add(nodeDistance);
        }
        pathLengthAccept.add(pathLength);
        nodeHeight0Accept.add(nodeHeights[0]);
        nodeHeight1Accept.add(nodeHeights[1]);
        leafCounts0Accept.add(leafCounts[0]);
        leafCounts1Accept.add(leafCounts[1]);
    }

    public void reject() {
        super.reject();

        if (mode != NARROW) {
            nodeDistanceReject.add(nodeDistance);
        }
        pathLengthReject.add(pathLength);
        nodeHeight0Reject.add(nodeHeights[0]);
        nodeHeight1Reject.add(nodeHeights[1]);
        leafCounts0Reject.add(leafCounts[0]);
        leafCounts1Reject.add(leafCounts[1]);
    }

    public LogColumn[] getColumns() {
        List<LogColumn> columns = new ArrayList<LogColumn>(Arrays.asList(super.getColumns()));

        if (mode != NARROW) {
            columns.add(getOperatorColumnInt("nodeDistAcc", nodeDistanceAccept));
            columns.add(getOperatorColumnInt("nodeDistRej", nodeDistanceReject));
        }
        columns.add(getOperatorColumnDouble("pathLengthAcc", pathLengthAccept));
        columns.add(getOperatorColumnDouble("pathLengthRej", pathLengthReject));
        columns.add(getOperatorColumnDouble("nodeheight0Acc", nodeHeight0Accept));
        columns.add(getOperatorColumnDouble("nodeheight0Rej", nodeHeight0Reject));
        columns.add(getOperatorColumnDouble("nodeheight1Acc", nodeHeight1Accept));
        columns.add(getOperatorColumnDouble("nodeheight1Rej", nodeHeight1Reject));
        columns.add(getOperatorColumnInt("leafCounts0Acc", leafCounts0Accept));
        columns.add(getOperatorColumnInt("leafCounts0Rej", leafCounts0Reject));
        columns.add(getOperatorColumnInt("leafCounts1Acc", leafCounts1Accept));
        columns.add(getOperatorColumnInt("leafCounts1Rej", leafCounts1Reject));

        return columns.toArray(new LogColumn[columns.size()]);
    }

    public double getMinimumAcceptanceLevel() {
        if( mode == NARROW ) {
            return 0.05;
        } else {
            return 0.01;
        }
    }

    public double getMinimumGoodAcceptanceLevel() {
        if( mode == NARROW ) {
            return 0.05;
        } else {
            return 0.01;
        }
    }

    public String getPerformanceSuggestion() {
        return "";

//        if( MCMCOperator.Utils.getAcceptanceProbability(this) < getMinimumAcceptanceLevel() ) {
//            return "";
//        } else if( MCMCOperator.Utils.getAcceptanceProbability(this) > getMaximumAcceptanceLevel() ) {
//            return "";
//        } else {
//            return "";
//        }
    }

//    public static final String INTERMEDIATE_EXCHANGE = "intermediateExchange";
/* The INTERMEDIATE_EXCHANGE is deprecated for unknown reason, therefore comment out its parser to make other 2 parsers registered properly
    public static XMLObjectParser INTERMEDIATE_EXCHANGE_PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return INTERMEDIATE_EXCHANGE;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            final TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
            final double weight = xo.getDoubleAttribute("weight");

            return new ExchangeOperator(INTERMEDIATE, treeModel, weight);
        }

        // ************************************************************************
        // AbstractXMLObjectParser implementation
        // ************************************************************************

        public String getParserDescription() {
            return "This element represents a intermediate exchange operator. "
                    + "This operator swaps two random subtrees.";
        }

        public Class getReturnType() {
            return ExchangeOperator.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newDoubleRule("weight"),
                new ElementRule(TreeModel.class)};
    }; */
}
