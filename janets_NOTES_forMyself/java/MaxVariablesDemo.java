/** 
 * MaxVariablesDemo.java is an application that compiles and runs
 * under J2SE 5.0. It requires no other files.
 */

public class MaxVariablesDemo {
    public static void main(String args[]) {

        //integers
        byte largestByte = Byte.MIN_VALUE;
        short largestShort = Short.MIN_VALUE;
        int largestInteger = Integer.MIN_VALUE;
        long largestLong = Long.MIN_VALUE;

        //real numbers
        float largestFloat = Float.MIN_VALUE;
        double largestDouble = Double.MIN_VALUE;

        //other primitive types
        char aChar = 'S';
        boolean aBoolean = false;

        //Display them all.
        System.out.println("The largest byte value is "
                           + largestByte + ".");
        System.out.println("The largest short value is "
                           + largestShort + ".");
        System.out.println("The largest integer value is "
                           + largestInteger + ".");
        System.out.println("The largest long value is "
                           + largestLong + ".");

        System.out.println("The largest float value is "
                           + largestFloat + ".");
        System.out.println("The largest double value is "
                           + largestDouble + ".");

        if (Character.isUpperCase(aChar)) {
            System.out.println("The character " + aChar
                               + " is uppercase.");
        } else {
            System.out.println("The character " + aChar
                               + " is lowercase.");
        }
        System.out.println("The value of aBoolean is "
                           + aBoolean + ".");
    }
}
