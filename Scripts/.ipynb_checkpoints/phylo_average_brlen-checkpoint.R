## File: phylo_average_brlen.R

## Usage: subtree.mean.branch.length <- phylo.average.brlen (subtree)

## Parameters: This function takes a tree in phylo format (from the ape package). This is usually a subtree with a root edge. 
## Value: A double representing the phylogenetic average of branch lengths for the clade.

## Example: 

## test.tree <- read.tree(read.tree(text= "(((A:1, B:1):1, C:1):1, D:1);"))
## plot(test.tree)
## test.subtree <- drop.tip(test.tree, "D", root.edge = 1)
## plot(test.subtree, root.edge=T)
## av.blen <- phylo.average.brlen(test.subtree)
## av.blen # should be 2.5


require(ape)

phylo.average.brlen <- function (phy)
{

    po <- reorder.phylo(phy, order = "postorder")

    br <- array(0, dim = Nnode(phy) + Ntip(phy))

    for (i in 1:Nedge(phy))
    {
        br[po$edge[i, 1]] <- br[po$edge[i, 1]] + (po$edge.length[i] + br[po$edge[i, 2]])/2
    }

    root.edge.length = ifelse(is.null(phy$root.edge), 0, phy$root.edge)
    singleton = ifelse(Ntip(phy) == 1, 2, 1) 

    return (br[Ntip(phy) + 1] * singleton + root.edge.length)

}



