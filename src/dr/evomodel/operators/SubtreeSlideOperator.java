/*
 * SubtreeSlideOperator.java
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
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.operators.SubtreeSlideOperatorParser;
import dr.inference.model.Statistic;
import dr.inference.operators.*;
import dr.math.MathUtils;

import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Implements the subtree slide move.
 *
 * @author Alexei Drummond
 * @version $Id: SubtreeSlideOperator.java,v 1.15 2005/06/14 10:40:34 rambaut Exp $
 */
public class SubtreeSlideOperator extends AbstractAdaptableTreeOperator {

    private static final boolean DEBUG = false;

    private DefaultTreeModel tree = null;
    private double size = 1.0;
    private boolean gaussian = false;
    private final boolean swapInRandomRate;
    private final boolean swapInRandomTrait;
    private final boolean scaledDirichletBranches;
    private AdaptationMode mode = AdaptationMode.DEFAULT;
    private final double targetAcceptance;

    private int nodeDistance;
    private double[] nodeHeights = new double[3];
    private int leafCounts;

    private List<Integer> nodeDistanceAccept;
    private List<Integer> nodeDistanceReject;
    private List<Double> nodeHeight0Accept;
    private List<Double> nodeHeight0Reject;
    private List<Double> nodeHeight1Accept;
    private List<Double> nodeHeight1Reject;
    private List<Double> nodeHeight2Accept;
    private List<Double> nodeHeight2Reject;
    private List<Integer> leafCountsAccept;
    private List<Integer> leafCountsReject;

    private int listMaxSize = 1000000;

    public SubtreeSlideOperator(DefaultTreeModel tree, double weight, double size, boolean gaussian,
                                boolean swapRates, boolean swapTraits, boolean scaleDirichletBranches,
                                AdaptationMode mode, double targetAcceptance) {
        super(mode, targetAcceptance);

        this.tree = tree;
        setWeight(weight);

        if (size == 0.0) {
            double b = 0.0;
            for (int k = 0; k < tree.getNodeCount(); ++k) {
                b += tree.getBranchLength(tree.getNode(k));
            }
            size = b / (2 * tree.getNodeCount());
        }

        this.size = size;
        this.gaussian = gaussian;
        this.swapInRandomRate = swapRates;
        this.swapInRandomTrait = swapTraits;
        this.scaledDirichletBranches = scaleDirichletBranches;

        this.mode = mode;
        this.targetAcceptance = targetAcceptance;

        nodeDistanceAccept = new ArrayList<Integer>();
        nodeDistanceReject = new ArrayList<Integer>();
        nodeHeight0Accept = new ArrayList<Double>();
        nodeHeight0Reject = new ArrayList<Double>();
        nodeHeight1Accept = new ArrayList<Double>();
        nodeHeight1Reject = new ArrayList<Double>();
        nodeHeight2Accept = new ArrayList<Double>();
        nodeHeight2Reject = new ArrayList<Double>();
        leafCountsAccept = new ArrayList<Integer>();
        leafCountsReject = new ArrayList<Integer>();
    }

    /**
     * Do a probablistic subtree slide move.
     *
     * @return the log-transformed hastings ratio
     */
    public double doOperation() {

        double logq;

        final NodeRef root = tree.getRoot();
        final double oldTreeHeight = tree.getNodeHeight(root);

        NodeRef i;

        // 1. choose a random node avoiding root
        do {
            i = tree.getNode(MathUtils.nextInt(tree.getNodeCount()));
        } while (root == i);

        final NodeRef iP = tree.getParent(i);
        final NodeRef CiP = getOtherChild(tree, iP, i);
        final NodeRef PiP = tree.getParent(iP);

        // 2. choose a delta to move
        final double delta = getDelta();
        final double oldHeight = tree.getNodeHeight(iP);
        final double newHeight = oldHeight + delta;

        nodeHeights[0] = tree.getNodeHeight(iP);
        nodeHeights[1] = newHeight;
        nodeHeights[2] = tree.getNodeHeight(i);
        leafCounts = TreeUtils.getLeafCount(tree, i);
        nodeDistance = 0;

        if (DEBUG) {
            System.out.println("\nSubTreeSlideOperator: oldTreeHeight = " + oldTreeHeight);
            System.out.println("Node selected: " + i + " ; delta = " + delta);
        }

        // 3. if the move is up
        if (delta > 0) {

            // 3.1 if the topology will change
            if (PiP != null && tree.getNodeHeight(PiP) < newHeight) {
                // find new parent
                NodeRef newParent = PiP;
                NodeRef newChild = iP;
                while (tree.getNodeHeight(newParent) < newHeight) {
                    newChild = newParent;
                    newParent = tree.getParent(newParent);
                    nodeDistance++;
                    if (newParent == null) break;
                }

                tree.beginTreeEdit();

                // 3.1.1 if creating a new root
                if (tree.isRoot(newChild)) {
                    tree.removeChild(iP, CiP);
                    tree.removeChild(PiP, iP);
                    tree.addChild(iP, newChild);
                    tree.addChild(PiP, CiP);
                    tree.setRoot(iP);
                    //System.err.println("Creating new root!");

                    if (tree.hasNodeTraits()) {
                        // **********************************************
                        // swap traits and rates so that root keeps it trait and rate values
                        // **********************************************

                        tree.swapAllTraits(newChild, iP);

                    }

                    if (tree.hasRates()) {
                        final double rootNodeRate = tree.getNodeRate(newChild);
                        tree.setNodeRate(newChild, tree.getNodeRate(iP));
                        tree.setNodeRate(iP, rootNodeRate);
                    }

                    // **********************************************

                }
                // 3.1.2 no new root
                else {
                    tree.removeChild(iP, CiP);
                    tree.removeChild(PiP, iP);
                    tree.removeChild(newParent, newChild);
                    tree.addChild(iP, newChild);
                    tree.addChild(PiP, CiP);
                    tree.addChild(newParent, iP);
                    //System.err.println("No new root!");
                }

                tree.setNodeHeight(iP, newHeight);

                tree.endTreeEdit();

                // 3.1.3 count the hypothetical sources of this destination.
                final int possibleSources = intersectingEdges(tree, newChild, oldHeight, null);
                //System.out.println("possible sources = " + possibleSources);

                logq = -Math.log(possibleSources);

            } else {
                // just change the node height
                tree.setNodeHeight(iP, newHeight);
                logq = 0.0;
            }
        }
        // 4 if we are sliding the subtree down.
        else {

            // 4.0 is it a valid move?
            if (tree.getNodeHeight(i) > newHeight) {
                nodeDistance = -1;
                return Double.NEGATIVE_INFINITY;
            }

            // 4.1 will the move change the topology
            if (tree.getNodeHeight(CiP) > newHeight) {

                List<NodeRef> newChildren = new ArrayList<NodeRef>();
                final int possibleDestinations = intersectingEdges(tree, CiP, newHeight, newChildren);

                if (DEBUG) {
                    System.out.println("possibleDestinations = " + possibleDestinations);
                }

                // if no valid destinations then return a failure
                if (newChildren.size() == 0) {
                    nodeDistance = -1;
                    return Double.NEGATIVE_INFINITY;
                }

                // pick a random parent/child destination edge uniformly from options
                final int childIndex = MathUtils.nextInt(newChildren.size());
                if (DEBUG) {
                    for (NodeRef ref : newChildren) {
                        System.out.println("child: " + ref.getNumber());
                    }
                }
                NodeRef newChild = newChildren.get(childIndex);
                NodeRef newParent = tree.getParent(newChild);
                nodeDistance = getNodeDistance(tree, iP, newParent);

                if (DEBUG) {
                    System.out.println("childIndex: " + childIndex);
                    System.out.println(newChild);
                    System.out.println(newParent);
                }

                tree.beginTreeEdit();

                // 4.1.1 if iP was root
                if (tree.isRoot(iP)) {
                    if (DEBUG) {
                        System.out.println("isRoot");
                    }
                    // new root is CiP
                    tree.removeChild(iP, CiP);
                    tree.removeChild(newParent, newChild);
                    tree.addChild(iP, newChild);
                    tree.addChild(newParent, iP);
                    tree.setRoot(CiP);

                    if (tree.hasNodeTraits()) {
                        // **********************************************
                        // swap traits and rates, so that root keeps it trait and rate values
                        // **********************************************

                        tree.swapAllTraits(iP, CiP);

                    }

                    if (tree.hasRates()) {
                        final double rootNodeRate = tree.getNodeRate(iP);
                        tree.setNodeRate(iP, tree.getNodeRate(CiP));
                        tree.setNodeRate(CiP, rootNodeRate);
                    }

                    // **********************************************

                    //System.err.println("DOWN: Creating new root!");
                } else {
                    tree.removeChild(iP, CiP);
                    tree.removeChild(PiP, iP);
                    tree.removeChild(newParent, newChild);
                    tree.addChild(iP, newChild);
                    tree.addChild(PiP, CiP);
                    tree.addChild(newParent, iP);
                    //System.err.println("DOWN: no new root!");
                }

                tree.setNodeHeight(iP, newHeight);

                tree.endTreeEdit();

                logq = Math.log(possibleDestinations);
            } else {
                tree.setNodeHeight(iP, newHeight);
                logq = 0.0;
            }
        }

        if (swapInRandomRate) {
            final NodeRef j = tree.getNode(MathUtils.nextInt(tree.getNodeCount()));
            if (j != i) {
                final double tmp = tree.getNodeRate(i);
                tree.setNodeRate(i, tree.getNodeRate(j));
                tree.setNodeRate(j, tmp);
            }

        }

        if (swapInRandomTrait) {
            final NodeRef j = tree.getNode(MathUtils.nextInt(tree.getNodeCount()));
            if (j != i) {

                tree.swapAllTraits(i, j);

//                final double tmp = tree.getNodeTrait(i, TRAIT);
//                tree.setNodeTrait(i, TRAIT, tree.getNodeTrait(j, TRAIT));
//                tree.setNodeTrait(j, TRAIT, tmp);
            }

        }

        // just return -Inf
        //if (logq == Double.NEGATIVE_INFINITY) throw new OperatorFailedException("invalid slide");

        if (scaledDirichletBranches) {
            if (oldTreeHeight != tree.getNodeHeight(tree.getRoot()))
                throw new UnsupportedOperationException("Temporarily disabled."); // TODO calculate Hastings ratio
        }

        return logq;
    }

    private double getDelta() {
        if (DEBUG) {
            System.out.println("size = " + size);
        }
        if (!gaussian) {
            return (MathUtils.nextDouble() * size) - (size / 2.0);
        } else {
            return MathUtils.nextGaussian() * size;
        }
    }

    private int intersectingEdges(Tree tree, NodeRef node, double height, List<NodeRef> directChildren) {

        final NodeRef parent = tree.getParent(node);

        /*if (DEBUG) {
            System.out.println("intersectingEdges");
            System.out.println("parent: " + parent);
        }*/

        if (tree.getNodeHeight(parent) < height) return 0;

        if (tree.getNodeHeight(node) < height) {
            if (directChildren != null) directChildren.add(node);
            return 1;
        }

        int count = 0;
        for (int i = 0; i < tree.getChildCount(node); i++) {
            count += intersectingEdges(tree, tree.getChild(node, i), height, directChildren);
        }
        return count;
    }

    public double getSize() {
        return size;
    }

    public void setSize(double size) {
        this.size = size;
    }

    @Override
    protected double getAdaptableParameterValue() {
        return Math.log(getSize());
    }

    @Override
    protected void setAdaptableParameterValue(double value) {
        setSize(Math.exp(value));
    }

    @Override
    public double getRawParameter() {
        return getSize();
    }

    @Override
    public String getAdaptableParameterName() {
        return "size";
    }

    public String getOperatorName() {
        return SubtreeSlideOperatorParser.SUBTREE_SLIDE + "(" + tree.getId() + ")";
    }

    public void accept(double deviation) {
        super.accept(deviation);

        nodeDistanceAccept.add(nodeDistance);
        nodeHeight0Accept.add(nodeHeights[0]);
        nodeHeight1Accept.add(nodeHeights[1]);
        nodeHeight2Accept.add(nodeHeights[2]);
        leafCountsAccept.add(leafCounts);

        if (nodeDistanceAccept.size() > listMaxSize) {
            nodeDistanceAccept.remove(0);
            nodeHeight0Accept.remove(0);
            nodeHeight1Accept.remove(0);
            nodeHeight2Accept.remove(0);
            leafCountsAccept.remove(0);
        }
    }

    public void reject() {
        super.reject();

        nodeDistanceReject.add(nodeDistance);
        nodeHeight0Reject.add(nodeHeights[0]);
        nodeHeight1Reject.add(nodeHeights[1]);
        nodeHeight2Reject.add(nodeHeights[2]);
        leafCountsReject.add(leafCounts);

        if (nodeDistanceReject.size() > listMaxSize) {
            nodeDistanceReject.remove(0);
            nodeHeight0Reject.remove(0);
            nodeHeight1Reject.remove(0);
            nodeHeight2Reject.remove(0);
            leafCountsReject.remove(0);
        }
    }

    public LogColumn[] getColumns() {
        List<LogColumn> columns = new ArrayList<LogColumn>(Arrays.asList(super.getColumns()));

        columns.add(getOperatorColumnInt("nodeDistAcc", nodeDistanceAccept));
        columns.add(getOperatorColumnInt("nodeDistRej", nodeDistanceReject));
        columns.add(getOperatorColumnDouble("nodeheightP0Acc", nodeHeight0Accept));
        columns.add(getOperatorColumnDouble("nodeheightP0Rej", nodeHeight0Reject));
        columns.add(getOperatorColumnDouble("nodeheightP1Acc", nodeHeight1Accept));
        columns.add(getOperatorColumnDouble("nodeheightP1Rej", nodeHeight1Reject));
        columns.add(getOperatorColumnDouble("nodeheightCAcc", nodeHeight2Accept));
        columns.add(getOperatorColumnDouble("nodeheightCRej", nodeHeight2Reject));
        columns.add(getOperatorColumnInt("leafCountsAcc", leafCountsAccept));
        columns.add(getOperatorColumnInt("leafCountsRej", leafCountsReject));

        return columns.toArray(new LogColumn[columns.size()]);
    }

}
