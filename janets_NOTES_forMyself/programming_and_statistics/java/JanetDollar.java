/** 
 * JanetDollar.java is an application that compiles and runs
 * under J2SE 5.0. It requires no other files.
 */

import java.text.DecimalFormat;

public class JanetDollar {
    public static void main(String[] args) {

        //a number
        float euronumber = 37;
        float conversion = 0.781162f;
        
        System.out.println("    euros = " + euronumber);

        System.out.println("converting...");
       
        float dollarnumber = euronumber / conversion  ;
 
        System.out.println("    dollars = " + dollarnumber);
        
        String pattern = "$###,###.##";
        DecimalFormat myFormatter = new DecimalFormat(pattern);
        String output = myFormatter.format(dollarnumber);

        System.out.println("    formatted dollars = " + output);
       
    }
}
