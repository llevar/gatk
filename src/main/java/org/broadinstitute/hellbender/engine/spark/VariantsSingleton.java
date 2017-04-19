package org.broadinstitute.hellbender.engine.spark;

import com.google.common.collect.Lists;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.utils.collections.IntervalsSkipList;
import org.broadinstitute.hellbender.utils.variant.GATKVariant;
import org.broadinstitute.hellbender.utils.variant.VariantContextVariantAdapter;

import java.io.File;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

public class VariantsSingleton {

    static String path;

    private final IntervalsSkipList<GATKVariant> intervalsSkipList;

    private VariantsSingleton() {
        System.out.println("Creating VariantsSingleton for " + path);
        // load from HDFS using NIO
        // ./src/test/resources/org/broadinstitute/hellbender/tools/BQSR/dbsnp_138.b37.excluding_sites_after_129.ch20.1m-1m1k.vcf
        List<GATKVariant> variants = getSerialVariants(path);
        intervalsSkipList = new IntervalsSkipList<>(variants);
        System.out.println("Finished creating VariantsSingleton for " + path);
    }

    /**
     * Loads variants using FeatureDataSource<VariantContext>.
     * @param vcf file to load
     * @return List of GATKVariants.
     */
    static List<GATKVariant> getSerialVariants(final String vcf) {
        try ( final FeatureDataSource<VariantContext> dataSource = new FeatureDataSource<>(new File(vcf), null, 0) ) {
            return Lists.newArrayList(wrapQueryResults(dataSource.iterator()));
        }
    }

    static List<GATKVariant> wrapQueryResults(final Iterator<VariantContext> queryResults ) {
        final List<GATKVariant> wrappedResults = new ArrayList<>();
        while ( queryResults.hasNext() ) {
            wrappedResults.add(VariantContextVariantAdapter.sparkVariantAdapter(queryResults.next()));
        }
        return wrappedResults;
    }

    // This class is only loaded when getInstance() is called.
    // JVM ensures it is only loaded once.
    // Each executor only has one copy, so it is shared by tasks.
    private static class SingletonHelper {
        private static final VariantsSingleton INSTANCE = new VariantsSingleton();
    }

    public static VariantsSingleton getInstance(String variantsPath) {
        path = variantsPath; // only effective the first time
        return SingletonHelper.INSTANCE;
    }

    public IntervalsSkipList<GATKVariant> getVariantsIntervalsSkipList() {
        return intervalsSkipList;
    }
}
