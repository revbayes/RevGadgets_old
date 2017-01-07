require(colorspace)

# modified from https://github.com/GuangchuangYu/ggtree/blob/master/R/tree-utilities.R
getXcoord2 <- function(x, root, parent, child, len, start=0, rev=FALSE) {
    x[root] <- start
    x[-root] <- NA  ## only root is set to start, by default 0

    currentNode <- root
    direction <- 1
    if (rev == TRUE) {
        direction <- -1
    }
    while(anyNA(x)) {
        idx <- which(parent %in% currentNode)
        newNode <- child[idx]
        x[newNode] <- x[parent[idx]]+len[idx] * direction
        currentNode <- newNode
    }

    return(x)
}

# modified from https://github.com/GuangchuangYu/ggtree/blob/master/R/tree-utilities.R
getXcoord <- function(tr) {
    edge <- tr$edge
    parent <- edge[,1]
    child <- edge[,2]
    root <- getRoot(tr)

    len <- tr$edge.length

    N <- getNodeNum(tr)
    x <- numeric(N)
    x <- getXcoord2(x, root, parent, child, len)
    return(x)
}

# modified from https://github.com/GuangchuangYu/ggtree/blob/master/R/tree-utilities.R
getYcoord <- function(tr, step=1) {
    Ntip <- length(tr[["tip.label"]])
    N <- getNodeNum(tr)

    edge <- tr[["edge"]]
    parent <- edge[,1]
    child <- edge[,2]

    cl <- split(child, parent)
    child_list <- list()
    child_list[as.numeric(names(cl))] <- cl

    y <- numeric(N)
    tip.idx <- child[child <= Ntip]
    y[tip.idx] <- 1:Ntip * step
    y[-tip.idx] <- NA

    currentNode <- 1:Ntip
    while(anyNA(y)) {
        pNode <- unique(parent[child %in% currentNode])
        ## piping of magrittr is slower than nested function call.
        ## pipeR is fastest, may consider to use pipeR
        ##
        ## child %in% currentNode %>% which %>% parent[.] %>% unique
        ## idx <- sapply(pNode, function(i) all(child[parent == i] %in% currentNode))
        idx <- sapply(pNode, function(i) all(child_list[[i]] %in% currentNode))
        newNode <- pNode[idx]

        y[newNode] <- sapply(newNode, function(i) {
            mean(y[child_list[[i]]], na.rm=TRUE)
            ##child[parent == i] %>% y[.] %>% mean(na.rm=TRUE)
        })

        currentNode <- c(currentNode[!currentNode %in% unlist(child_list[newNode])], newNode)
        ## currentNode <- c(currentNode[!currentNode %in% child[parent %in% newNode]], newNode)
        ## parent %in% newNode %>% child[.] %>%
        ##     `%in%`(currentNode, .) %>% `!` %>%
        ##         currentNode[.] %>% c(., newNode)
    }

    return(y)
}

# modified from https://github.com/GuangchuangYu/ggtree/blob/master/R/tree-utilities.R
getParent <- function(tr, node) {
    if ( node == getRoot(tr) )
        return(0)
    edge <- tr[["edge"]]
    parent <- edge[,1]
    child <- edge[,2]
    res <- parent[child == node]
    if (length(res) == 0) {
        stop("cannot found parent node...")
    }
    if (length(res) > 1) {
        stop("multiple parent found...")
    }
    return(res)
}

assign_state_labels = function(t, state_labels, include_start_states, n_states=1)
{

    # exit if no state labels provided
    if (is.null(state_labels)) {
        return(t)
    }
        

    # what is the ancestral state name tag?
    if (include_start_states) {
        state_pos_str_base = c("start_state_", "end_state_")
    } else {
        state_pos_str_base = c("anc_state_")
    }
    
    # create list of ancestral state name tags
    state_pos_str_to_update = c(sapply(1:n_states, function(x) { paste(state_pos_str_base,x,sep="")}))
    
    # overwrite state labels
    for (m in state_pos_str_to_update)
    {
        x_state = attributes(t)$stats[[m]]
        #print(x_state)
        #print(as.numeric(levels(x_state)))
        to_states = state_labels[as.numeric(levels(x_state))]
        if (any(is.na(to_states)))
        {
            #na_idx = which(is.na(to_states))
            #max_state = max(as.vector(to_states), na.rm=T)
            #print(na_idx)
            #print(max_state)
            #to_states[na_idx] = max_state+1:length(na_idx)
            #print(to_states)
            to_states = state_labels #[1:length(state_labels)]
        }
        
        x_state = plyr:::mapvalues(x_state, from=levels(x_state), to=to_states)
        x_state <- factor(x_state, levels=state_labels, ordered=T)
        #print(x_state)
        attributes(t)$stats[[m]] = x_state
    }
    return(t)
}

build_state_probs = function(t, state_labels, include_start_states) {

    n_states = length(state_labels)
    n_tips = length(attributes(t)$phylo$tip.label)
    n_node = 2 * n_tips - 1
    
    dat = list()
    
    if (include_start_states) {
        state_tags = c("start","end")
    } else {
        state_tags = c("anc")
    }
    
    for (s in state_tags) {
        dat[[s]] = data.frame( matrix(0, nrow=n_node, ncol=n_states) )
        dat[[s]] = cbind(node=1:n_node, dat[[s]])
        #colnames(dat[[s]]) = c("node", state_labels)
        
        for (i in 1:3)
        {
            m = paste(s,"_state_",i,sep="")
            pp_str = paste(m,"_pp",sep="")
            n_tmp = attributes(t)$stats$node
            #x_tmp = as.numeric(as.vector(attributes(t)$stats[[m]]))
            x_tmp = as.vector(attributes(t)$stats[[m]])
            pp_tmp = as.vector(attributes(t)$stats[[pp_str]])
            
            #print(m)
            ##print(pp_str)
            #print(attributes(t)$stats)
            ##print(n_tmp)
            #print(x_tmp)
            #print(pp_tmp)
            for (j in 1:length(x_tmp))
            {
                #cat(c(n_tmp[j],x_tmp[j],pp_tmp[j]),"\n",sep="\t")
                if (!is.na(x_tmp[j])) {
                    k = which(x_tmp[j]==state_labels)
                    dat[[s]][n_tmp[j], k+1] = pp_tmp[j]
                    #cat(c(n_tmp[j],x_tmp[j],pp_tmp[j]),"\n",sep="\t")
                }
            }
        }
        #print(head(dat[[s]]))
        #print(apply(dat[[s]][2:ncol(dat_start)], 1, sum))
    
    }
    
    
    return(dat)
}


################################################################################
#
# @brief Function to plot ancestral states and the associated uncertainty
#        for continuous and discrete characters.
#
#        For discrete characters 'summary_statistic="MAP"' should be used,
#        and for continuous characters 'summary_statistic="mean"'. If the
#        tree and tip labels do not fit in the screen, adjust the visible
#        area using the 'xlim_visible' argument.
#
#        If 'summary_statistic="MAP"', the maximum a posteriori ancestral
#        state will be plotted on each node. The color corresponds to the
#        character state, and the size of the circle represents the posterior
#        probability of that state. Cladogenetic models that estimate 
#        ancestral states for both the beginning and end of each branch
#        are plotted by setting "include_start_states=TRUE".
#
#        Maximum a posteriori ancestral chromosome numbers can be plotted 
#        with 'summary_statistic="MAPChromosome"'. For chromosomes the
#        color represents the posterior probability and the size of the
#        circle represents the chromosome number.
#
#        For 'summary_statistic="mean"' the color represents the size of 
#        the 95% confidence interval, and the size of the cirlce represents
#        represents the mean character state.
#
# @date Last modified: 2016-09-29
# @author Will Freyman
# @version 1.0
# @since 2016-08-31, version 1.0.0
#
# @param    tree_file               character     Path to the ancestral state tree generated by RevBayes.
# @param    summary_statistic       character   The type of summary statistic to plot.
# @param    tree_layout             character   One of 'rectangular', 'slanted', 'fan', 'circular', 'radial', or 'unrooted'.
# @param    include_start_states    logical     Plot start and end ancestral states. Used for cladogenetic models.
#
#
################################################################################
plot_ancestral_states = function(tree_file, 
                                 summary_statistic="MAP", 
                                 tree_layout="rectangular",
                                 include_start_states=FALSE, 
                                 xlim_visible=c(0, 40), 
                                 ylim_visible=NULL,
                                 tip_label_size=4, 
                                 tip_label_offset=3,
                                 tip_label_italics=FALSE,
                                 tip_node_size=2,
                                 node_label_size=4, 
                                 node_label_nudge_x=0.1, 
                                 shoulder_label_size=3, 
                                 shoulder_label_nudge_x=-0.1, 
                                 alpha=0.5, 
                                 node_size_range=c(6, 15), 
                                 color_low="#D55E00",
                                 color_mid="#F0E442",
                                 color_high="#009E73",
                                 show_state_legend=TRUE,
                                 show_posterior_legend=TRUE,
                                 show_tree_scale=TRUE,
                                 state_labels=NULL,
                                 state_colors=NULL,
                                 ...) { 

    if ( (summary_statistic %in% c("MAP", "mean", "MAPChromosome", "MAPRange", "PieRange")) == FALSE ) {
        print("Invalid summary statistic.")
        return()
    }

    # read in tree
    t = read.beast(tree_file)
    
    # add state labels
    t = assign_state_labels(t, state_labels, include_start_states)
    
    # add state colors
    use_state_colors = !is.null(state_colors)
    if (!is.null(state_colors) && !is.null(state_labels))
    {
        names(state_colors) = state_labels
    }
    
    # get anc state matrices (for pie/bar charts)
    dat_state = build_state_probs(t, state_labels, include_start_states)
    
    tree = attributes(t)$phylo
    n_node = getNodeNum(tree)

    # remove underscores from tip labels
    attributes(t)$phylo$tip.label = gsub("_", " ", attributes(t)$phylo$tip.label)
    
    if (tip_label_italics) {
        attributes(t)$phylo$tip.label = paste("italic('", attributes(t)$phylo$tip.label, "')", sep="")
    }

    # add tip labels
    p = ggtree(t, layout=tree_layout) 
    p = p + geom_tiplab(size=tip_label_size, offset=tip_label_offset, parse=tip_label_italics)
       
    if (summary_statistic == "MAPChromosome") {
        
        if (include_start_states) {
            
            if (!("start_state_1" %in% colnames(attributes(t)$stats))) {
                print("Start states not found in input tree.")
                return()
            }

            # add ancestral states as node labels
            p = p + geom_text(aes(label=end_state_1), hjust="left", nudge_x=node_label_nudge_x, size=node_label_size)

            # set the root's start state to NA
            attributes(t)$stats$start_state_1[n_node] = NA

            # add clado daughter lineage start states on "shoulders" of tree
            # get x, y coordinates of all nodes
            x = getXcoord(tree)
            y = getYcoord(tree)
            x_anc = numeric(n_node)
            node_index = numeric(n_node)
            for (i in 1:n_node) {
                if (getParent(tree, i) != 0) {
                    # if not the root, get the x coordinate for the parent node
                    x_anc[i] = x[getParent(tree, i)]
                    node_index[i] = i
                }
            }
            shoulder_data = data.frame(node=node_index, x_anc=x_anc, y=y)
            p = p %<+% shoulder_data
            
            # plot the states on the "shoulders"
            p = p + geom_text(aes(label=start_state_1, x=x_anc, y=y), hjust="right", nudge_x=shoulder_label_nudge_x, size=shoulder_label_size, na.rm=TRUE)
            
            # show ancestral states as size / posteriors as color
            p = p + geom_nodepoint(aes(colour=end_state_1_pp, size=end_state_1), alpha=alpha)
            min_low = 0.0
            max_up = 1.0
            p = p + scale_colour_gradient2(low=color_low, mid=color_mid, high=color_high, limits=c(min_low, max_up), midpoint=0.5)
            if (show_state_legend) {
                p = p + guides(size=guide_legend("Chromosome Number"))
            } else {
                p = p + guides(size=FALSE)
            }
            if (show_posterior_legend) {
                p = p + guides(colour=guide_legend("Posterior Probability", override.aes = list(size=8)))
            } else {
                p = p + guides(colour=FALSE)
            }
        }
    } else if (summary_statistic == "MAPRange") {
        if (!include_start_states) {
            warning("Ignoring that include_start_states is set to FALSE")
        }
        if (!("start_state_1" %in% colnames(attributes(t)$stats))) {
            print("Start states not found in input tree.")
            return()
        }

        # add ancestral states as node labels
        p = p + geom_text(aes(label=end_state_1), hjust="left", nudge_x=node_label_nudge_x, size=node_label_size)

        # set the root's start state to NA
        attributes(t)$stats$start_state_1[n_node] = NA

        # add clado daughter lineage start states on "shoulders" of tree
        # get x, y coordinates of all nodes
        x = getXcoord(tree)
        y = getYcoord(tree)
        x_anc = numeric(n_node)
        node_index = numeric(n_node)
        for (i in 1:n_node) {
            if (getParent(tree, i) != 0) {
                # if not the root, get the x coordinate for the parent node
                x_anc[i] = x[getParent(tree, i)]
                node_index[i] = i
            }
        }
        shoulder_data = data.frame(node=node_index, x_anc=x_anc, y=y)
        p = p %<+% shoulder_data
       
        # plot the states on the "shoulders"
        p = p + geom_text(aes(label=start_state_1, x=x_anc, y=y), hjust="right", nudge_x=shoulder_label_nudge_x, size=shoulder_label_size, na.rm=TRUE)
        p = p + geom_nodepoint(aes(colour=factor(start_state_1), x=x_anc, y=y, size=end_state_1_pp),na.rm=TRUE, alpha=alpha)
        p = p + geom_tippoint(aes(colour=factor(start_state_1), x=x_anc, y=y, size=end_state_1_pp),na.rm=TRUE, alpha=alpha)
        
        # show tip states as color
        #print(shoulder_data)
        #print(x_anc)
        #print(c(attributes(t)$stats$start_state_1,attributes(t)$stats$end_state_1))
       
        p = p + geom_tippoint(aes(colour=factor(end_state_1)), size=tip_node_size, alpha=alpha) 
        
        # show ancestral states as color / posteriors as size
        p = p + geom_nodepoint(aes(colour=factor(end_state_1), size=end_state_1_pp), alpha=alpha)

        if (show_state_legend) {
            p = p + guides(colour=guide_legend("Range", override.aes = list(size=8), order=1))
        } else {
            p = p + guides(colour=FALSE)
        }
        
        if (show_posterior_legend) {
            p = p + guides(size=guide_legend("Posterior probability", order=2))
        } else {
            p = p + guides(size=FALSE)
        }

    } else if (summary_statistic == "PieRange") {
        if (!include_start_states) {
            warning("Ignoring that include_start_states is set to FALSE")
        }
        if (!("start_state_1" %in% colnames(attributes(t)$stats))) {
            print("Start states not found in input tree.")
            return()
        }

        # add ancestral states as node labels
        p = p + geom_text(aes(label=end_state_1), hjust="left", nudge_x=node_label_nudge_x, size=node_label_size)

        # set the root's start state to NA
        attributes(t)$stats$start_state_1[n_node] = NA
        
        # add clado daughter lineage start states on "shoulders" of tree
        # get x, y coordinates of all nodes
        x = getXcoord(tree)
        y = getYcoord(tree)
        x_anc = numeric(n_node)
        node_index = numeric(n_node)
        for (i in 1:n_node) {
            if (getParent(tree, i) != 0) {
                # if not the root, get the x coordinate for the parent node
                x_anc[i] = x[getParent(tree, i)]
                node_index[i] = i
            }
        }
        shoulder_data = data.frame(node=node_index, x_anc=x_anc, y=y)
        p = p %<+% shoulder_data
       
        # plot the states on the "shoulders"
        #p = p + geom_text(aes(label=start_state_1, x=x_anc, y=y), hjust="right", nudge_x=shoulder_label_nudge_x, size=shoulder_label_size, na.rm=TRUE)
        
        # show ancestral states as color / posteriors as size
        #p = p + geom_nodepoint(aes(colour=factor(end_state_1), size=end_state_1_pp), alpha=alpha)
        
        pies = nodebar(dat_state[["start"]], cols=2:ncol(dat_state[["start"]]))
        inset(p,pies,width=1,height=1)
        
        # show tip states as color
        p = p + geom_tippoint(aes(colour=factor(end_state_1)), size=tip_node_size, alpha=alpha)
        
        p = p + theme_tree2()
        
        if (show_state_legend) {
            p = p + guides(size=guide_legend("Posterior probability"))
        } else {
            p = p + guides(size=FALSE)
        }
        if (show_posterior_legend) {
            p = p + guides(colour=guide_legend("Range", override.aes = list(size=8)))
        } else {
            p = p + guides(colour=FALSE)
        }

    } else if (summary_statistic == "MAP") {

        if (include_start_states) {
            print("Start states not yet implemented for MAP ancestral states.")
            return()
    
        }
        if (!("anc_state_1" %in% colnames(attributes(t)$stats))) {
            anc_data = data.frame(node=names(attributes(t)$stats$end_state_1), 
                                  anc_state_1=levels(attributes(t)$stats$end_state_1)[attributes(t)$stats$end_state_1],
                                  anc_state_1_pp=as.numeric(levels(attributes(t)$stats$end_state_1_pp))[attributes(t)$stats$end_state_1_pp])
            p = p %<+% anc_data
        }
        
        # add ancestral states as node labels
        p = p + geom_text(aes(label=anc_state_1), hjust="left", nudge_x=node_label_nudge_x, size=node_label_size)

        # show ancestral states as color / posteriors as size
        p = p + geom_nodepoint(aes(colour=factor(anc_state_1), size=anc_state_1_pp), alpha=alpha)
        
        # show the tip values
        p = p + geom_tippoint(aes(colour=factor(anc_state_1)), size=tip_node_size, alpha=alpha)
        
        # set up the legend
        if (show_state_legend) {
            p = p + guides(colour=guide_legend("State"))        
        } else {
            p = p + guides(colour=FALSE)
        }
        if (show_posterior_legend) {
            p = p + guides(size=guide_legend("Posterior Probability"))
        } else {
            p = p + guides(size=FALSE)
        }

    } else if (summary_statistic == "mean") {
    
        if (include_start_states) {
            print("Start states not implemented for mean ancestral states.")
            return()
        }

        # add ancestral states as node labels
        p = p + geom_text(aes(label=round(mean, 2)), hjust="left", nudge_x=node_label_nudge_x, size=node_label_size)

        # show the size of the 95% CI as color 
        lowers = as.numeric(levels(attributes(t)$stats$lower_0.95_CI))[attributes(t)$stats$lower_0.95_CI]
        uppers = as.numeric(levels(attributes(t)$stats$upper_0.95_CI))[attributes(t)$stats$upper_0.95_CI]
        diffs = uppers - lowers
        diffs_df = data.frame(node=names(attributes(t)$stats$lower_0.95_CI), diff_vals=diffs)
        p = p %<+% diffs_df 

        min_low = min(diffs, na.rm=TRUE)
        max_up = max(diffs, na.rm=TRUE)
        mid_val = min_low + (max_up - min_low) / 2.0
        p = p + scale_colour_gradient2(low=color_low, mid=color_mid, high=color_high, limits=c(min_low, max_up), midpoint=mid_val)
        p = p + geom_nodepoint(aes(size=mean, colour=diff_vals), alpha=alpha)

        # show the tip values
        p = p + geom_tippoint(aes(size=mean), color="grey", alpha=alpha)

        # set up the legend
        if (show_state_legend) {
            legend_text = "Mean State"
            p = p + guides(size=guide_legend(legend_text))
        } else {
            p = p + guides(size=FALSE)
        }
        if (show_posterior_legend) {
            p = p + guides(colour=guide_legend("95% CI Width", override.aes=list(size=4)))
        } else {
            p = p + guides(colour=FALSE)
        }
    
    } 
    
    
    if (use_state_colors) {
        p = p + scale_color_manual(values=state_colors, breaks=state_labels)
    }
    
    p = p + scale_radius(range = node_size_range)
    p = p + theme(legend.position="left")
    
    if (show_tree_scale)
    {
       #p = p + theme_tree2()
    }

    # set visible area
    p = p + coord_cartesian(xlim = xlim_visible, ylim=ylim_visible, expand=TRUE)
    print(p)
    return(p)
}

