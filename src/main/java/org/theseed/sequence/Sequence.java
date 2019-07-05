/**
 *
 */
package org.theseed.sequence;

/**
 * This class represents a record from a FASTA file.  It contains a sequence,
 * a label, and a comment.
 *
 * @author Bruce Parrello
 *
 */
public class Sequence {

    // FIELDS
    String label;
    String comment;
    String sequence;

    /**
     * Construct a blank, empty sequence.
     */
    public Sequence() {
        this.label = "";
        this.comment = "";
        this.sequence = "";
    }


    /**
     * Construct a sequence from known strings.
     *
     * @param label		the label to assign
     * @param comment	the command about the sequence
     * @param sequence	the sequence text
     */
    public Sequence(String label, String comment, String sequence) {
        super();
        this.label = label;
        this.comment = comment;
        this.sequence = sequence;
    }

    /**
     * @return the length of the sequence
     */
    public int length() {
        return this.sequence.length();
    }

    /**
     * @return the label
     */
    public String getLabel() {
        return label;
    }

    /**
     * @param label the label to set
     */
    public void setLabel(String label) {
        this.label = label;
    }

    /**
     * @return the comment
     */
    public String getComment() {
        return comment;
    }

    /**
     * @param comment the comment to set
     */
    public void setComment(String comment) {
        this.comment = comment;
    }

    /**
     * @return the sequence
     */
    public String getSequence() {
        return sequence;
    }

    /**
     * @param sequence the sequence to set
     */
    public void setSequence(String sequence) {
        this.sequence = sequence;
    }

}
