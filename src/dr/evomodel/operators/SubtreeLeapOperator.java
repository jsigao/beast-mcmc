/*
 * SubtreeLeapOperator.java
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

package dr.evomodel.operators;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.operators.SubtreeLeapOperatorParser;
import dr.evomodelxml.operators.TipLeapOperatorParser;
import dr.inference.distribution.CauchyDistribution;
import dr.inference.loggers.LogColumn;
import dr.inference.operators.AdaptationMode;
import dr.math.MathUtils;
import dr.math.distributions.Distribution;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * Implements the Subtree Leap move.
 *
 * This move picks a node at random (except for the root) and then moves the parent to any location
 * that is a certain patristic distance from its starting point (the distance is drawn from a Gaussian).
 *
 * It is always possible for the node to move up (potentially becoming the root) but the destination can't
 * be younger than the original node. All possible destinations are collected and then picked amongst
 * uniformly.
 *
 * @author Andrew Rambaut
 * @author Luiz Max Carvalho
 * @author Mathieu Fourment
 */
public class SubtreeLeapOperator extends AbstractAdaptableTreeOperator {

    public enum DistanceKernelType {
        NORMAL("normal") {
            @Override
            double getDelta(double size) {
                return Math.abs(MathUtils.nextGaussian() * size);
            }
        },
        CAUCHY("Cauchy") {
            @Override
            double getDelta(double size) {
                Distribution distK = new CauchyDistribution(0, size);
                double u = MathUtils.nextDouble();
                return Math.abs(distK.quantile(u));
            }
        };

        DistanceKernelType(String name) {
            this.name = name;
        }

        @Override
        public String toString() {
            return name;
        }

        String name;

        abstract double getDelta(double size);
    }

    private double size;

    private final TreeModel tree;
    private final DistanceKernelType distanceKernel;
    private final boolean slideOnly;

    private final List<Integer> tips;

    private int nodeDistance;
    private double pathLength;
    private double[] nodeHeights = new double[2];
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

    /**
     * Constructor
     *
     * @param tree   the tree
     * @param weight the weight
     * @param size   scaling on a unit Gaussian to draw the patristic distance from
     * @param targetAcceptance the desired acceptance probability
     * @param distanceKernel the distribution from which to draw the patristic distance
     * @param mode   coercion mode
     */
    public SubtreeLeapOperator(TreeModel tree, double weight, double size, DistanceKernelType distanceKernel, AdaptationMode mode, double targetAcceptance) {
        this(tree, weight, size, distanceKernel, false, mode, targetAcceptance);
    }

    /**
     * Constructor
     *
     * @param tree   the tree
     * @param weight the weight
     * @param size   scaling on a unit Gaussian to draw the patristic distance from
     * @param targetAcceptance the desired acceptance probability
     * @param distanceKernel the distribution from which to draw the patristic distance
     * @param slideOnly if true, only slide up and down the tree, never across (mimics SubtreeSlide)
     * @param mode   coercion mode
     */
    public SubtreeLeapOperator(TreeModel tree, double weight, double size, DistanceKernelType distanceKernel, boolean slideOnly, AdaptationMode mode, double targetAcceptance) {
        super(mode, targetAcceptance);

        this.tree = tree;
        setWeight(weight);
        this.size = size;
        this.distanceKernel = distanceKernel;
        this.slideOnly = slideOnly;
        this.tips = null;

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

    /**
     * Constructor that takes a taxon set to pick from for the move.
     *
     * @param tree   the tree
     * @param taxa   some taxa
     * @param weight the weight
     * @param size   scaling on a unit Gaussian to draw the patristic distance from
     * @param mode   coercion mode
     */
    public SubtreeLeapOperator(TreeModel tree, TaxonList taxa, double weight, double size, DistanceKernelType distanceKernel, AdaptationMode mode, double targetAcceptance) {
        super(mode, targetAcceptance);

        this.tree = tree;
        setWeight(weight);
        this.size = size;
        this.distanceKernel = distanceKernel;
        this.slideOnly = false;
        this.tips = new ArrayList<Integer>();

        for (Taxon taxon : taxa) {
            boolean found = false;
            for (int i = 0; i < tree.getExternalNodeCount(); i++) {
                NodeRef tip = tree.getExternalNode(i);
                if (tree.getNodeTaxon(tip).equals(taxon)) {
                    tips.add(tip.getNumber());
                    found = true;
                    break;
                }
            }

            if (!found) {
                throw new IllegalArgumentException("Taxon, " + taxon.getId() + ", not found in tree with id " + tree.getId());
            }
        }

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


    /**
     * Do a subtree leap move.
     *
     * @return the log-transformed hastings ratio
     */
    public double doOperation() {
        double logq;

        final double delta = distanceKernel.getDelta(size);

        final NodeRef root = tree.getRoot();

        NodeRef node;

        if (tips == null) {
            // Pick a node (but not the root)
            do {
                // choose a random node avoiding root
                node = tree.getNode(MathUtils.nextInt(tree.getNodeCount()));

            } while (node == root);
        } else {
            // Pick a tip from the specified set of tips.
            node = tree.getNode(tips.get(MathUtils.nextInt(tips.size())));
        }

        // get its parent - this is the node we will prune/graft
        final NodeRef parent = tree.getParent(node);

        // get the node's sibling
        final NodeRef sibling = getOtherChild(tree, parent, node);

        // and its grand parent
        final NodeRef grandParent = tree.getParent(parent);

        final Map<NodeRef, Double> destinations = getDestinations(node, parent, sibling, delta, slideOnly);
        final List<NodeRef> destinationNodes = new ArrayList<NodeRef>(destinations.keySet());

        // pick uniformly from this list
        int r = MathUtils.nextInt(destinations.size());

        double forwardProbability = 1.0 / destinations.size();

        final NodeRef j = destinationNodes.get(r);
        final double newHeight = destinations.get(j);

        final NodeRef jParent = tree.getParent(j);

        if (logOperatorStat) {
            NodeRef ca = TreeUtils.getCommonAncestorNode(tree, parent, j);
            if (ca == j) {
                nodeDistance = getNodeDistance(tree, parent, j);
            } else {
                nodeDistance = getNodeDistance(tree, parent, j) - 1;
            }
            pathLength = delta;
            nodeHeights[0] = tree.getNodeHeight(parent);
            nodeHeights[1] = newHeight;
            if (logCladeOperated) {
                cladeIndices[0] = getCladeIdx(tree, node);
                cladeIndices[1] = getCladeIdx(tree, j);
            }
        }

        if (jParent != null && newHeight > tree.getNodeHeight(jParent)) {
            throw new IllegalArgumentException("height error");
        }

        if (newHeight < tree.getNodeHeight(j)) {
            throw new IllegalArgumentException("height error");
        }

        tree.beginTreeEdit();

        if (j == parent || jParent == parent) {
            // the subtree is not actually moving but the height will change
        } else {
            if (grandParent == null) {
                // if the parent of the original node is the root then the sibling becomes
                // the root.
                tree.removeChild(parent, sibling);
                tree.setRoot(sibling);

            } else {
                // remove the parent of node by connecting its sibling to its grandparent.
                tree.removeChild(parent, sibling);
                tree.removeChild(grandParent, parent);
                tree.addChild(grandParent, sibling);
            }

            if (jParent == null) {
                // adding the node to the root of the tree
                tree.addChild(parent, j);
                tree.setRoot(parent);
            } else {
                // remove destination edge j from its parent
                tree.removeChild(jParent, j);

                // add destination edge to the parent of node
                tree.addChild(parent, j);

                // and add the parent of i as a child of the former parent of j
                tree.addChild(jParent, parent);
            }
        }

        tree.setNodeHeight(parent, newHeight);

        tree.endTreeEdit();

        if (tree.getParent(parent) != null && newHeight > tree.getNodeHeight(tree.getParent(parent))) {
            throw new IllegalArgumentException("height error");
        }

        if (newHeight < tree.getNodeHeight(node)) {
            throw new IllegalArgumentException("height error");
        }

        if (newHeight < tree.getNodeHeight(getOtherChild(tree, parent, node))) {
            throw new IllegalArgumentException("height error");
        }

        final Map<NodeRef, Double> reverseDestinations = getDestinations(node, parent, getOtherChild(tree, parent, node), delta, slideOnly);
        double reverseProbability = 1.0 / reverseDestinations.size();

        // hastings ratio = reverse Prob / forward Prob
        logq = Math.log(reverseProbability) - Math.log(forwardProbability);
        return logq;
    }

    private Map<NodeRef, Double> getDestinations(NodeRef node, NodeRef parent, NodeRef sibling, double delta, boolean slideOnly) {

        final Map<NodeRef, Double> destinations = new LinkedHashMap<NodeRef, Double>();

        // get the parent's height
        final double height = tree.getNodeHeight(parent);

        final double heightBelow = height - delta;

        if (heightBelow > tree.getNodeHeight(node)) {
            // the destination height below the parent is compatible with the node
            // see if there are any destinations on the sibling's branch
            final List<NodeRef> edges = new ArrayList<NodeRef>();

            getIntersectingEdges(tree, sibling, heightBelow, edges);

            // add the intersecting edges and the height
            for (NodeRef n : edges) {
                destinations.put(n, heightBelow);
            }
        }

        final double heightAbove = height + delta;

        NodeRef node1 = parent;

        // walk up to root
        boolean done = false;
        while (!done) {
            NodeRef parent1 = tree.getParent(node1);

            if (parent1 != null) {
                final double height1 = tree.getNodeHeight(parent1);
                if (height1 < heightAbove) {
                    if (!slideOnly) { // if we are not just sliding up or down...
                        
                        // We haven't reached the height above the original height so go down
                        // the sibling subtree to look for other possible destinations
                        NodeRef sibling1 = getOtherChild(tree, parent1, node1);

                        double heightBelow1 = height1 - (heightAbove - height1);

                        if (heightBelow1 > tree.getNodeHeight(node)) {

                            final List<NodeRef> edges = new ArrayList<NodeRef>();

                            getIntersectingEdges(tree, sibling1, heightBelow1, edges);

                            // add the intersecting edges and the height
                            for (NodeRef n : edges) {
                                destinations.put(n, heightBelow1);
                            }
                        }
                    }
                } else {
                    // add the current node as a destination
                    destinations.put(node1, heightAbove);
                    done = true;
                }

                node1 = parent1;
            } else {
                // node1 is the root - add it as a destination and stop loop
                destinations.put(node1, heightAbove);
                done = true;
            }
        }

        return destinations;
    }

    private int getIntersectingEdges(Tree tree, NodeRef node, double height, List<NodeRef> edges) {

        final NodeRef parent = tree.getParent(node);

        if (tree.getNodeHeight(parent) < height) return 0;

        if (tree.getNodeHeight(node) < height) {
            edges.add(node);
            return 1;
        }

        int count = 0;
        for (int i = 0; i < tree.getChildCount(node); i++) {
            count += getIntersectingEdges(tree, tree.getChild(node, i), height, edges);
        }
        return count;
    }

    public double getSize() {
        return size;
    }

    public void setSize(double size) {
        assert Double.isFinite(size);
        this.size = size;
    }

    @Override
    protected void setAdaptableParameterValue(double value) {
        setSize(Math.exp(value));
    }

    @Override
    protected double getAdaptableParameterValue() {
        return Math.log(getSize());
    }

    @Override
    public double getRawParameter() {
        return getSize();
    }

    public String getAdaptableParameterName() {
        return "size";
    }

    public String getOperatorName() {
        if (tips == null) {
            return SubtreeLeapOperatorParser.SUBTREE_LEAP + "(" + tree.getId() + ")";
        } else {
            return TipLeapOperatorParser.TIP_LEAP + "(" + tree.getId() + ")";
        }
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
