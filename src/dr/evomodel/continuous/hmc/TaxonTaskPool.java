package dr.evomodel.continuous.hmc;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

/**
 * @author Marc A. Suchard
 * @author Andrew Holbrook
 */

public class TaxonTaskPool {

    class TaxonTaskIndices {

        final int start;
        final int stop;
        final int task;

        TaxonTaskIndices(int start, int stop, int task) {
            this.start = start;
            this.stop = stop;
            this.task = task;
        }

        public String toString() {
            return start + " " + stop;
        }
    }

    final private ExecutorService pool;
    final private List<TaxonTaskIndices> indices;

    TaxonTaskPool(int taxonCount, int threadCount) {
        this.indices = setupTasks(taxonCount, Math.abs(threadCount));
        this.pool = setupParallelServices(threadCount);
    }

//    public ExecutorService getPool() { return pool; }

//    public List<TaxonTaskIndices> getIndices() { return indices; }

    public int getNumThreads() { return indices.size(); }

    private List<TaxonTaskIndices> setupTasks(int taxonCount, int threadCount) {
        List<TaxonTaskIndices> tasks = new ArrayList<TaxonTaskIndices>(threadCount);

        int length = taxonCount / threadCount;
        if (taxonCount % threadCount != 0) ++length;

        int start = 0;

        for (int task = 0; task < threadCount && start < taxonCount; ++task) {
            tasks.add(new TaxonTaskIndices(start, Math.min(start + length, taxonCount), task));
            start += length;
        }

        return tasks;
    }

    private ExecutorService setupParallelServices(int threadCount) {

        final ExecutorService pool;

        if (threadCount > 1) {
            pool = Executors.newFixedThreadPool(threadCount);
        } else if (threadCount < 0) {
            pool = Executors.newCachedThreadPool();
        } else {
            pool = null;
        }

        return pool;
    }

    interface TaxonCallable {
        void execute(int taxon, int thread);
    }

    public void fork(final TaxonCallable runnable) {
        if (indices.size() == 1) {

            final TaxonTaskIndices index = indices.get(0);
            for (int taxon = index.start; taxon < index.stop; ++taxon) {
                runnable.execute(taxon, 0);
            }

        } else {

            List<Callable<Object>> calls = new ArrayList<Callable<Object>>();

            for (final TaxonTaskIndices indexSet : indices) {

                calls.add(Executors.callable(
                        new Runnable() {
                            @Override
                            public void run() {
                                for (int taxon = indexSet.start; taxon < indexSet.stop; ++taxon) {
                                    runnable.execute(taxon, indexSet.task);
                                }
                            }
                        }
                ));
            }

            try {
                pool.invokeAll(calls);
            } catch (InterruptedException exception) {
                exception.printStackTrace();
            }
        }
    }
}
