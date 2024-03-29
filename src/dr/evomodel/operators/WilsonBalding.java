/*
 * WilsonBalding.java
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
import dr.evomodelxml.operators.WilsonBaldingParser;
import dr.math.MathUtils;

import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Implements the unweighted wilson-balding branch swapping move.
 *
 * @author Alexei Drummond
 * @version $Id: WilsonBalding.java,v 1.38 2005/06/14 10:40:34 rambaut Exp $
 */
public class WilsonBalding extends AbstractTreeOperator implements ConstrainableTreeOperator {

    private TreeModel tree = null;

    private int nodeDistance;
    private double pathLength;
    private double[] nodeHeights = new double[2];
    private int leafCounts;
    private int[] cladeIndices = new int[2];
    protected boolean logOperatorStat = false;

    private List<Integer> nodeDistanceAccept;
    private List<Integer> nodeDistanceReject;
    private List<Double> pathLengthAccept;
    private List<Double> pathLengthReject;
    private List<Double> nodeHeight0Accept;
    private List<Double> nodeHeight0Reject;
    private List<Double> nodeHeight1Accept;
    private List<Double> nodeHeight1Reject;
    private List<Integer> cladeIndex0Accept;
    private List<Integer> cladeIndex0Reject;
    private List<Integer> cladeIndex1Accept;
    private List<Integer> cladeIndex1Reject;
    private List<Long> calculationCountAccept;
    private List<Long> calculationCountReject;

    public WilsonBalding(TreeModel tree, double weight) {
        this.tree = tree;
        setWeight(weight);

        nodeDistanceAccept = new ArrayList<Integer>();
        nodeDistanceReject = new ArrayList<Integer>();
        pathLengthAccept = new ArrayList<Double>();
        pathLengthReject = new ArrayList<Double>();
        nodeHeight0Accept = new ArrayList<Double>();
        nodeHeight0Reject = new ArrayList<Double>();
        nodeHeight1Accept = new ArrayList<Double>();
        nodeHeight1Reject = new ArrayList<Double>();
        cladeIndex0Accept = new ArrayList<Integer>();
        cladeIndex0Reject = new ArrayList<Integer>();
        cladeIndex1Accept = new ArrayList<Integer>();
        cladeIndex1Reject = new ArrayList<Integer>();
        calculationCountAccept = new ArrayList<Long>();
        calculationCountReject = new ArrayList<Long>();
    }

    public double doOperation(TreeModel tree) {

        double logq = proposeTree(tree);
        int tipCount = tree.getExternalNodeCount();
        if (tree.getExternalNodeCount() != tipCount) {
            int newCount = tree.getExternalNodeCount();
            throw new RuntimeException("Lost some tips in modified SPR! (" +
                    tipCount + "-> " + newCount + ")");
        }
        //System.out.println("last accepted deviation: " + getDeviation());
        //System.out.println("logq=" + logq);
        return logq;
    }

    public double doOperation(){
        return doOperation(tree);
    }

    /**
     * WARNING: Assumes strictly bifurcating tree.
     */
    public double proposeTree(TreeModel tree){

        NodeRef i;
        double oldMinAge, newMinAge, newRange, oldRange, newAge, q;

        //Bchoose

        //for (int n =0; n < tree.getNodeCount(); n++) {
        //	System.out.println(n + " " + ( (tree.getNode(n) == null) ? "null" : tree.getNode(n).getId()));
        //}

        // choose a random node avoiding root
        final int nodeCount = tree.getNodeCount();
        do {
            i = tree.getNode(MathUtils.nextInt(nodeCount));
        } while (tree.getRoot() == i);
        final NodeRef iP = tree.getParent(i);

        // choose another random node to insert i above
        NodeRef j = tree.getNode(MathUtils.nextInt(nodeCount));
        NodeRef k = tree.getParent(j);

        // make sure that the target branch <k, j> is above the subtree being moved
        while ((k != null && tree.getNodeHeight(k) <= tree.getNodeHeight(i)) || (i == j)) {
            j = tree.getNode(MathUtils.nextInt(nodeCount));
            k = tree.getParent(j);
        }

        NodeRef ca = TreeUtils.getCommonAncestorNode(tree, iP, j);
        if (logOperatorStat) {
            if (ca == j) {
                nodeDistance = getNodeDistance(tree, iP, j);
            } else {
                nodeDistance = getNodeDistance(tree, iP, j) - 1;
            }
            pathLength = -1;
            nodeHeights[0] = tree.getNodeHeight(iP);
            nodeHeights[1] = -1;
            if (logCladeOperated) {
                cladeIndices[0] = getCladeIdx(tree, i);
                cladeIndices[1] = getCladeIdx(tree, j);
            }
        }

        // disallow moves that change the root.
        if (j == tree.getRoot() || iP == tree.getRoot()) {
            return Double.NEGATIVE_INFINITY;
        }

        if (k == iP || j == iP || k == i) {
            return Double.NEGATIVE_INFINITY;
        }

        final NodeRef CiP = getOtherChild(tree, iP, i);
        NodeRef PiP = tree.getParent(iP);

//		ConstantPopulation demoFunc = null;
//		if (demoModel != null && demoModel.getDemographicFunction() instanceof ConstantPopulation) {
//			demoFunc = (ConstantPopulation)demoModel.getDemographicFunction();
//		}

//        if (j == tree.getRoot()) {
//			if (demoModel != null) {
//				delta = -demoFunc.getN0() * Math.log(MathUtils.nextDouble());
//			} else {
//				delta = tree.getNodeHeight(j) * MathUtils.nextDouble();
//			}
//			newAge = tree.getNodeHeight(j) + delta;
//
//			PiP = tree.getParent(iP);
//			oldMinAge = Math.max(tree.getNodeHeight(i), tree.getNodeHeight(CiP));
//			oldRange = tree.getNodeHeight(PiP) - oldMinAge;
//
//			if (demoFunc == null) {
//				q = tree.getNodeHeight(j) / oldRange;
//			} else {
//				q = Math.exp(delta/demoFunc.getN0())*demoFunc.getN0()/oldRange;
//			}
//		} else if (iP == tree.getRoot()) {
//
//			newMinAge = Math.max(tree.getNodeHeight(i), tree.getNodeHeight(j));
//			newRange = tree.getNodeHeight(k) - newMinAge;
//			newAge = newMinAge + (MathUtils.nextDouble()*newRange);
//
//			if (demoFunc == null) {
//				if (tree.getNodeHeight(iP) > (tree.getNodeHeight(CiP) * 2)) throw new OperatorFailedException("too big");
//				q = newRange / tree.getNodeHeight(CiP);
//			} else {
//				q = (tree.getNodeHeight(CiP)-tree.getNodeHeight(iP))/demoFunc.getN0() + Math.log(newRange/demoFunc.getN0());
//			}
//		} else {
        newMinAge = Math.max(tree.getNodeHeight(i), tree.getNodeHeight(j));
        newRange = tree.getNodeHeight(k) - newMinAge;
        newAge = newMinAge + (MathUtils.nextDouble() * newRange);
        oldMinAge = Math.max(tree.getNodeHeight(i), tree.getNodeHeight(CiP));
        oldRange = tree.getNodeHeight(PiP) - oldMinAge;
        q = newRange / Math.abs(oldRange);

        if (logOperatorStat) {
            nodeHeights[1] = newAge;
            if (ca == j) {
                pathLength = newAge - tree.getNodeHeight(iP);
            } else {
                pathLength = TreeUtils.getPathLength(tree, iP, k) + tree.getNodeHeight(k) - newAge;
            }
        }
        
        //System.out.println(newRange + "/" + oldRange + "=" + q);
//		}

        //Bupdate

        tree.beginTreeEdit();

        if (j == tree.getRoot()) {

            // 1. remove edges <iP, CiP>
            tree.removeChild(iP, CiP);
            tree.removeChild(PiP, iP);

            // 2. add edges <k, iP>, <iP, j>, <PiP, CiP>
            tree.addChild(iP, j);
            tree.addChild(PiP, CiP);

            // iP is the new root
            tree.setRoot(iP);

        } else if (iP == tree.getRoot()) {

            // 1. remove edges <k, j>, <iP, CiP>, <PiP, iP>
            tree.removeChild(k, j);
            tree.removeChild(iP, CiP);

            // 2. add edges <k, iP>, <iP, j>, <PiP, CiP>
            tree.addChild(iP, j);
            tree.addChild(k, iP);

            //CiP is the new root
            tree.setRoot(CiP);

        } else {
            // 1. remove edges <k, j>, <iP, CiP>, <PiP, iP>
            tree.removeChild(k, j);
            tree.removeChild(iP, CiP);
            tree.removeChild(PiP, iP);

            // 2. add edges <k, iP>, <iP, j>, <PiP, CiP>
            tree.addChild(iP, j);
            tree.addChild(k, iP);
            tree.addChild(PiP, CiP);
        }

        tree.setNodeHeight(iP, newAge);

        tree.endTreeEdit();

        // AR - I don't believe this check is needed and in tests it never fails...
//        try {
//            tree.checkTreeIsValid();
//        } catch( MutableTree.InvalidTreeException ite ) {
//            throw new RuntimeException(ite.toString());
////            throw new OperatorFailedException(ite.toString());
//        }


        return Math.log(q);
    }

    public double getMinimumAcceptanceLevel() {
        return 0.01;
    }

    public String getPerformanceSuggestion() {
        // seems like equvivalent code to me :)
        return "";

//        if (MCMCOperator.Utils.getAcceptanceProbability(this) < getMinimumAcceptanceLevel()) {
//            return "";
//        } else if (MCMCOperator.Utils.getAcceptanceProbability(this) > getMaximumAcceptanceLevel()) {
//            return "";
//        } else {
//            return "";
//        }
    }

    public String getOperatorName() {
        return WilsonBaldingParser.WILSON_BALDING + "(" + tree.getId() + ")";
    }

    public void accept(double deviation) {
        super.accept(deviation);

        if (logOperatorStat) {
            nodeDistanceAccept.add(nodeDistance);
            pathLengthAccept.add(pathLength);
            nodeHeight0Accept.add(nodeHeights[0]);
            nodeHeight1Accept.add(nodeHeights[1]);
            if (logCladeOperated) {
                cladeIndex0Accept.add(cladeIndices[0]);
                cladeIndex1Accept.add(cladeIndices[1]);
            }
            calculationCountAccept.add(calculationCount);
        }
    }

    public void reject() {
        super.reject();

        if (logOperatorStat) {
            nodeDistanceReject.add(nodeDistance);
            pathLengthReject.add(pathLength);
            nodeHeight0Reject.add(nodeHeights[0]);
            nodeHeight1Reject.add(nodeHeights[1]);
            if (logCladeOperated) {
                cladeIndex0Reject.add(cladeIndices[0]);
                cladeIndex1Reject.add(cladeIndices[1]);
            }
            calculationCountReject.add(calculationCount);
        }
    }

    public LogColumn[] getColumns() {
        List<LogColumn> columns = new ArrayList<LogColumn>(Arrays.asList(super.getColumns()));
        logOperatorStat = true;

        columns.add(getOperatorColumnInt("nodeDistAcc", nodeDistanceAccept));
        columns.add(getOperatorColumnInt("nodeDistRej", nodeDistanceReject));
        columns.add(getOperatorColumnDouble("pathLengthAcc", pathLengthAccept));
        columns.add(getOperatorColumnDouble("pathLengthRej", pathLengthReject));
        columns.add(getOperatorColumnDouble("nodeheightP0Acc", nodeHeight0Accept));
        columns.add(getOperatorColumnDouble("nodeheightP0Rej", nodeHeight0Reject));
        columns.add(getOperatorColumnDouble("nodeheightP1Acc", nodeHeight1Accept));
        columns.add(getOperatorColumnDouble("nodeheightP1Rej", nodeHeight1Reject));

        if (logCladeOperated) {
            columns.add(getOperatorColumnInt("cladeId0Acc", cladeIndex0Accept));
            columns.add(getOperatorColumnInt("cladeId0Rej", cladeIndex0Reject));
            columns.add(getOperatorColumnInt("cladeId1Acc", cladeIndex1Accept));
            columns.add(getOperatorColumnInt("cladeId1Rej", cladeIndex1Reject));
        }

        columns.add(getOperatorColumnLong("calcountAcc", calculationCountAccept));
        columns.add(getOperatorColumnLong("calcountRej", calculationCountReject));

        return columns.toArray(new LogColumn[columns.size()]);
    }
}
